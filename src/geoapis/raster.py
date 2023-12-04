# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 10:10:55 2022

@author: pearsonra
"""
import abc
import logging
import pathlib
import time
import typing
import urllib
import zipfile

import geopandas
import geopandas as gpd
import requests
from tqdm import tqdm

from geoapis.vector import Linz as LinzVector


class KoordinatesExportsQueryBase(abc.ABC):
    """An abstract class to manage fetching Raster data using the Koordinates exports
    API. Downloads the GeoTiff specified in the run routine.

    API details at: https://help.koordinates.com/site-admin-apis/export-api/

    Parameters
    ----------

    key: str
        The API key.  Must have Search, view and download exported data
        permissions.

    cache_path: pathlib.Path
        The location to download all GeoTiffs queried in the run method.

    crs: int
        The CRS EPSG code for the GeoTifss to be downloaded as.

    bounding_polygon: geopandas.geodataframe.GeoDataFrame
        An option geometry to clip the downloaded GeoTiffs within.

    """

    @property
    @abc.abstractmethod
    def NETLOC_API():
        """This should be instantiated in the base class. Provide the netloc of the data
        service."""

        raise NotImplementedError("NETLOC_API must be instantiated in the child class")

    SCHEME = "https"
    PATH = "services/api/v1"
    PATH_API_END = "/exports"
    K_CRS = "EPSG:4326"

    def __init__(
        self,
        key: str,
        cache_path: typing.Union[str, pathlib.Path],
        crs: int = None,
        bounding_polygon: geopandas.GeoSeries = None,
    ):
        """Load in the wfs key and CRS/bounding_polygon if specified. Specify the layer
        to import during run."""

        self.key = key
        self.cache_path = pathlib.Path(cache_path)
        self.bounding_polygon = bounding_polygon
        self.crs = crs

        self.base_url = urllib.parse.urlunparse(
            (
                self.SCHEME,
                self.NETLOC_API,
                self.PATH,
                "",
                "",
                "",
            )
        )

        self._set_up()

    def _set_up(self):
        """Ensure the bouding_polygon and CRS are in agreement."""

        # Set the crs from the bounding_polygon if it's not been set
        if self.crs is None and self.bounding_polygon is not None:
            logging.info("The download CRS is being set from the bounding_polygon")
            self.crs = self.bounding_polygon.crs.to_epsg()
        # Set the bounding_polygon crs from the crs if they differ
        if (
            self.bounding_polygon is not None
            and self.crs != self.bounding_polygon.crs.to_epsg()
        ):
            logging.info(
                "The bounding_polygon is being transformed to the specified "
                "download CRS"
            )
            self.bounding_polygon.to_crs(self.crs)
        # Enforce the bounding_polygon must be a single geometry if it exists
        if self.bounding_polygon is not None:
            self.bounding_polygon = self.bounding_polygon.explode(index_parts=False)
            if not (self.bounding_polygon.type == "Polygon").all():
                logging.warning(
                    "All bounding_polygon parts aren't Polygon's. Ignoring"
                    f" those that aren't {self.bounding_polygon.geometry}"
                )

                self.bounding_polygon = self.bounding_polygon[
                    self.bounding_polygon.type == "Polygon"
                ]
            number_of_coords = sum(
                [
                    len(polygon.coords)
                    for polygon in self.bounding_polygon.explode(
                        index_parts=False
                    ).exterior
                ]
            )
            assert number_of_coords < 1000, (
                "The bounding polygon must be less than 1000 points. Consider using the"
                " bbox to simplify the geometry"
            )

    def run(self, layer: int, index_tiles: int = None) -> None:
        """Query for a specified layer and return a geopandas.GeoDataFrame of the vector
        features. If a polygon_boundary is specified, only return vectors passing
        through this polygon."""

        headers = {"Authorization": f"key {self.key}"}

        # Create the initial request
        api_query = {
            "crs": f"EPSG:{self.crs}",
            "formats": {"grid": "image/tiff;subtype=geotiff"},
            "items": [{"item": f"{self.base_url}/layers/{layer}/"}],
        }
        if self.bounding_polygon is not None:
            polygon_chuck = (
                gpd.GeoSeries(self.bounding_polygon.unary_union.buffer(0), crs=2193)
                .explode(index_parts=False)
                .reset_index(drop=True)
            )
            polygon_chuck = polygon_chuck.to_crs(self.K_CRS)
            if (
                len(polygon_chuck) == 1
            ):  # member must be an array of LinearRing coordinate arrays
                exterior = polygon_chuck.loc[0].exterior
                query_coords = [list(exterior.coords)]
                api_query["extent"] = {
                    "type": "Polygon",
                    "coordinates": query_coords,
                }
        logging.info("Send initial request to download image")
        response = requests.post(
            url=f"{self.base_url}/exports/", headers=headers, json=api_query
        )
        try:
            query_id = response.json()["id"]
            name = response.json()["name"]
            self.export(layer, name, query_id, headers)
        except KeyError:
            logging.warning("The query failed. Check the invalid_reasons")
            vec = LinzVector(
                key=self.key,
                crs=self.crs,
                bounding_polygon=self.bounding_polygon,
            )
            index_tiles = vec.get_features(index_tiles)
            index_tiles.sort_values(by="tilename", inplace=True)
            # split into chucks of 10000
            n = 5000  # chunk row size
            if len(index_tiles) < n:
                n = int(len(index_tiles) / 2) + 1
            index_chucks = [
                index_tiles[i : i + n] for i in range(0, index_tiles.shape[0], n)
            ]
            for index_chunk in index_chucks:
                polygon_chuck = (
                    gpd.GeoSeries(index_chunk.unary_union.buffer(0), crs=2193)
                    .explode(index_parts=False)
                    .reset_index(drop=True)
                )
                polygon_chuck = polygon_chuck.to_crs(self.K_CRS)
                if (
                    len(polygon_chuck) == 1
                ):  # member must be an array of LinearRing coordinate arrays
                    exterior = polygon_chuck.loc[0].exterior
                    query_coords = [list(exterior.coords)]
                    api_query["extent"] = {
                        "type": "Polygon",
                        "coordinates": query_coords,
                    }
                else:  # member must be an array of Polygon coordinate arrays
                    query_coords = [
                        [list(polygon.exterior.coords)] for polygon in polygon_chuck
                    ]
                    api_query["extent"] = {
                        "type": "MultiPolygon",
                        "coordinates": query_coords,
                    }
                response = requests.post(
                    url=f"{self.base_url}/exports/", headers=headers, json=api_query
                )
                query_id = response.json()["id"]
                name = response.json()["name"]
                self.export(layer, name, query_id, headers)
            pass

    def export(self, layer: int, name: str, query_id: int, headers: dict) -> None:
        # Check the state of your exports until the triggered raster exports completes
        logging.info("Check status of download request")
        with tqdm(
            total=1.0,
            desc="Generating export...",
            bar_format="{l_bar}{bar}| {elapsed}<{remaining}",
        ) as pbar:
            pbar.update(0)
            while True:
                response = requests.get(
                    f"{self.base_url}/exports/",
                    headers=headers,
                )
                # find the triggered export
                element = [
                    element for element in response.json() if element["id"] == query_id
                ][0]
                logging.info(f"/texport state is {element['state']}")
                if element["state"] == "processing":
                    progress_url = element["url"]
                    with requests.get(
                        progress_url,
                        headers=headers,
                    ) as progress_response:
                        progress_response.raise_for_status()
                        progress = progress_response.json()["progress"]
                        if progress:
                            pbar.update(progress - pbar.n)
                        else:
                            pbar.update(1 - pbar.n)
                    logging.info("Not complete - check again in 20s")
                    time.sleep(10)
                    continue
                elif element["state"] == "complete":
                    pbar.update(1 - pbar.n)
                    logging.info("/tCompleted - move to download")
                    break
                else:
                    logging.warning(
                        f"Could not download raster. Ended with status {element['state']}"
                    )
                    return
        # Download the completed export
        logging.info(f"Downloading {element['download_url']} to {self.cache_path}")
        zip_path = self.cache_path / f"{layer}.zip"
        with open(zip_path, mode="wb") as zip_file:
            with requests.get(
                element["download_url"],
                headers={"Authorization": f"key {self.key}"},
                stream=True,
            ) as response:
                response.raise_for_status()
                total = int(response.headers.get("content-length", 0))
                tqdm_params = {
                    "desc": "Downloading raster",
                    "total": total,
                    "miniters": 1,
                    "unit": "B",
                    "unit_scale": True,
                    "unit_divisor": 1024,
                }
                with tqdm(**tqdm_params) as progress_bar:
                    for chunk in response.iter_content(chunk_size=1024):
                        progress_bar.update(len(chunk))
                        zip_file.write(chunk)

        with zipfile.ZipFile(zip_path) as zip_object:
            # extract with progress bar
            tqdm_params = {
                "desc": "Extracting raster",
                "total": len(zip_object.infolist()),
                "miniters": 1,
            }
            if name.startswith("lds-"):
                name = name[4:]
            if name.endswith("-GTiff"):
                name = name[:-6]
            for file in tqdm(zip_object.infolist(), **tqdm_params):
                # check if the file exists
                if (self.cache_path / f"{name}" / file.filename).exists():
                    logging.warning(
                        f"File {self.cache_path / f'{name}' / file.filename} already exists"
                    )
                    continue
                zip_object.extract(file, self.cache_path / f"{name}")


class Linz(KoordinatesExportsQueryBase):
    """A class to manage fetching Vector data from LINZ.

    LIRS data service can be accessed at: https://https://data.linz.govt.nz/

    Note that only rasters supporting the grid image/tiff geotiff are supported
    """

    NETLOC_API = "data.linz.govt.nz"


class Lris(KoordinatesExportsQueryBase):
    """A class to manage fetching Vector data from LRIS.

    LIRS data service can be accessed at: https://lris.scinfo.org.nz/

    Note that only rasters supporting the grid image/tiff geotiff are supported
    """

    NETLOC_API = "lris.scinfo.org.nz"


class StatsNz(KoordinatesExportsQueryBase):
    """A class to manage fetching Vector data from the Stats NZ datafinder.

    Stats NZ data service can be accessed at: datafinder.stats.govt.nz

    Note that only rasters supporting the grid image/tiff geotiff are supported
    """

    NETLOC_API = "datafinder.stats.govt.nz"


class KoordinatesQuery(KoordinatesExportsQueryBase):
    """A class to manage fetching Vector data from any generic data portal supporting
    WFS.

    Note that the 'geometry_name' used when making a WFS 'cql_filter' queries can vary
    between layers. You will need to specify the 'geometry_name' of the layers you want
    to download.
    """

    def __init__(
        self,
        key: str,
        netloc_url: str,
        crs: int = None,
        bounding_polygon: geopandas.geodataframe.GeoDataFrame = None,
    ):
        """Set NETLOC_API and instantiate the KoordinatesExportsQueryBase"""

        self.netloc_url = netloc_url

        # Setup the WfsQueryBase class
        super(KoordinatesQuery, self).__init__(
            key=key, crs=crs, bounding_polygon=bounding_polygon
        )

    @property
    def NETLOC_API(self):
        """Instantiate the entered netloc of the data service."""

        return self.netloc_url
