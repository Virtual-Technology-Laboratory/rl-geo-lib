# Copyright (c) 2015, Roger Lew (rogerlew.gmail.com)
# Date: 12/9/2014
# License: BSD (3-clause license)
# 
# The project described was supported by NSF award number IIA-1301792
# from the NSF Idaho EPSCoR Program and by the National Science Foundation.

from __future__ import print_function

import concurrent.futures
from concurrent.futures import ThreadPoolExecutor

import os
import warnings
import sys
import urllib

from ast import literal_eval
from collections import namedtuple
from csv import DictReader
from glob import glob
from Queue import Queue

from os.path import join as joinpath
from zipfile import ZipFile

class Tile(object):
    def __init__(self, name, extents, url):
        extents = literal_eval(extents)
        xmin = min(extents[0], extents[2])
        xmax = max(extents[0], extents[2])
        ymin = min(extents[1], extents[3])
        ymax = max(extents[1], extents[3])
        self.extents = [xmin, ymin, xmax, ymax]
        
        self.name = name
        self.url = url

    def intersects(self, other_extents):
        o_xmin = min(other_extents[0], other_extents[2])
        o_xmax = max(other_extents[0], other_extents[2])
        o_ymin = min(other_extents[1], other_extents[3])
        o_ymax = max(other_extents[1], other_extents[3])
        
        s_xmin, s_ymin, s_xmax, s_ymax = self.extents

        if o_xmax <= s_xmin:
            return False
        if o_xmin >= s_xmax:
            return False
        if o_ymax <= s_ymin:
            return False
        if o_ymin >= s_ymax:
            return False

        return True
         
def _downloadTile(tupledArgs):
    imgDB, name, overwrite, verbose = tupledArgs
    
    db_location = imgDB.db_location
    dest_dir = imgDB._extractPath(name)
    datasetFN = imgDB._datasetFn(name)

    if os.path.exists(datasetFN) and not overwrite:
        warnings.warn('"%s" has already been downloaded' % dest_dir)
        return
    
    tile = imgDB.manifest[name]
    url = tile.url
    zip_fn = joinpath(db_location, name + '.zip')

    if verbose:
        print('Retrieving "%s"' % url, end='')

    def reporthook(count, blockSize, totalSize):
        # insideIdaho is returning  a totalSize of -1
        if count % 1000 == 0:
            sys.stdout.write(".") # every 8 MB
            sys.stdout.flush()

    (filename, headers) = urllib.urlretrieve(url, zip_fn, reporthook)
    
    if verbose:
        print('\nExtracting "%s"' % zip_fn)
    with ZipFile(zip_fn) as zf:
        zf.extractall(dest_dir)

    if verbose:
        print('Removing "%s"\n' % zip_fn)
    os.remove(zip_fn)

    return 1

        
class ImageDB(object):
    def __init__(self, db_location, manifest_fn):
        if not os.path.exists(db_location):
            os.mkdir(db_location)
            
        self.db_location = db_location
        self.manifest_fn = manifest_fn

        # set self.manifest
        self._readManifest()

    def _readManifest(self):
        manifest = {}
        
        f = open(self.manifest_fn)
        rdr = DictReader(f)
        for L in rdr:
            manifest[L['name']] = Tile(**L)
        f.close()

        self.manifest = manifest

    def _getAcquired(self):
        wc = joinpath(self.db_location, '*')
        acquired = []

        for fn in glob(wc):
            if os.path.isdir(fn):
                name = os.path.basename(fn)
                if os.path.exists(self._datasetFn(name)):
                    acquired.append(name)
            
        return acquired

    def _getNeeded(self, extents):
        needed = []
        for tile in self.manifest.values():
            if tile.intersects(extents):
                needed.append(tile.name)
        return needed
                
    def fetchForExtents(self, extents, overwrite=False,
                        verbose=True, simultaneous=4):
        """
        Syncronizes the database to the extents and return a list
        of paths to datasets intersecting the extents

        Parameters
        ----------
        extents : listable of floats
            Can be ordered [xmin, ymin, xmax, ymax] or
            [xmin, ymax, xmax, ymin]

        verbose : boolean
            Specifies whether to print feedback during the process

        Returns
        -------
        dslist : list
            list containing paths to the dataset files
        """
        needed = self._getNeeded(extents)
        n = len(needed)
            
        if n == 0:
            raise Exception('Extents are outside the database')

        if verbose:
            print('Need %i tile%s.' % (n, ('s', '')[n == 1]))
        
        acquired = self._getAcquired()
        queue = [fn for fn in needed if fn not in acquired]
        q = len(queue)
        
        if verbose:
            print('Need to acquire %i tiles.' % q)

        tupledArgs = []
        for name in queue:
            tupledArgs.append(
                (self, name, overwrite, verbose))

        futures = Queue()
        with ThreadPoolExecutor(max_workers=simultaneous) as executor:

            for args in tupledArgs:
                future = executor.submit(_downloadTile, args)

##        while not futures.empty():
##            future = futures.get()
##            print('future',future.done())
            
##        for tup in tupledArgs:
##            _downloadTile(tup)

        return self.getDatasetList(extents)
        
    def getDatasetList(self, extents=None):
        """
        Returns a list to the paths of all the acquired
        datasets
        """
        if extents is None:
            acquired = self._getAcquired()        
        else:
            acquired = self._getNeeded(extents)        

        dslist = [self._datasetFn(name) for name in acquired]
        
        for fn in dslist:
            if not os.path.exists(fn):
                raise Exception('Cannot locate "%s"' % fn)

        return dslist

class OrthoImagery1m(ImageDB):
    """
    Connector for the insideIdaho.org 1m Orthoimagery data
    """
    def __init__(self, db_location):
        dirname = os.path.dirname(__file__)
        self.db_manifest = joinpath(dirname, '1mOrthoimagery.db.csv')
        
        super(OrthoImagery1m, self)\
                .__init__(db_location, self.db_manifest)

    def _datasetFn(self, name):
        """
        Specifies the location of the raster dataset relative
        to the db_location.
        """    
        return joinpath(self.db_location, name, name + '.tif')

    def _extractPath(self, name):
        """
        Specifies the directory to extract downloaded datasets
        relative to the db_location
        """    
        return joinpath(self.db_location, name)


class NED10m(ImageDB):
    """
    Connector for the insideIdaho.org 10m USGS National Elevation
    data.
    """
    def __init__(self, db_location):
        dirname = os.path.dirname(__file__)
        self.db_manifest = joinpath(dirname, '10mNED.db.csv')
        
        super(NED10m, self)\
                .__init__(db_location, self.db_manifest)

    def _datasetFn(self, name):
        """
        Specifies the absolute path to the location of the raster
        dataset relative.
        """    
        return joinpath(self.db_location, name, name, 'w001001.adf')

    def _extractPath(self, name):
        """
        Specifies the directory to extract downloaded datasets
        relative to the db_location
        """    
        return self.db_location


_url_ortho_template = \
    'http://cloud.insideidaho.org/arcgis/rest/services/'\
    'imageryBaseMapsEarthCover/2011_1m_idaho/ImageServer/'\
    'exportImage?bbox={bbox}&bboxSR={bboxSR}&size={size}'\
    '&imageSR={imageSR}&time=&format=tiff&pixelType=U8&noData=0'\
    '&noDataInterpretation=esriNoDataMatchAny'\
    '&interpolation=+RSP_BilinearInterpolation'\
    '&compression=&compressionQuality=&bandIds=&mosaicRule='\
    '&renderingRule=&f=pjson'

def fetch_orthoimagery(extents, dst_fname,
                       bboxSR=4296, size=(4096,4096), imageSR=4296):
    """
    fetches orthoimagery from insideidaho.org REST Image Export Service

    extents should be ordered xmin, ymin, xmax, ymax
    """
    
    # limited to predefined wkids. Could not get this to work with custom
    # projection.
    
    global _url_ortho_template
    
    bbox = '%2C'.join(map(str, extents))
    bboxSR = str(bboxSR)     # bounding box spatial reference
    size = '%i%%2C%i' % size # requested size in pixels
    imageSR = str(imageSR)   # spatial reference of requested image

    url = _url_ortho_template.format(bbox=bbox, bboxSR=bboxSR,
                                     size=size, imageSR=imageSR)

    print('  Requesting JSON...')
    r = requests.get(url)
    if r.status_code != 200:
        warnings.warn('url request returned %i' % r.status_code)
        return

    d = literal_eval(r.text)

    print('  Fetching Image...')
    urllib.urlretrieve (d['href'], dst_fname)
    
if __name__ == '__main__':
    
    extents = '[-116.4633,43.8147,-116.1227,44.0111]'
    extents = literal_eval(extents)

    db = OrthoImagery1m(r'D:\tmp')    
    db.fetchForExtents(extents)

##    NED_db = NED10m(r'D:\ownCloud\documents\geo_data\idaho_ned10m')
##    
##    NED_db.fetchForExtents(extents)
##    pprint(NED_db.getDatasetList())
