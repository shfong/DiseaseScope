import requests

class BigGim(object):
    """Client for BigGIM
    
    Inspired by this repository: 
    https://github.com/NCATS-Tangerine/BigGIM_APIWrapper/blob/master/app.py
    that wraps the BigGIM API.
    """

    URL = "http://biggim.ncats.io/api"

    def __init__(self):
        pass

    def _getBG(self, endpoint, data=None):
        """GET method for BigGIM"""

        if data is None: 
            data = {}

        try: 
            req = requests.get(f"{self.URL}/{endpoint}", data=data)
            req.raise_for_status()

            return req.json()

        except requests.HTTPError as e: 
            return {'error': str(e)}

    def _postBG(self, endpoint, data=None): 
        """POST method for BigGIM"""

        if data is None: 
            data = {}

        req = requests.post(f"{self.URL}/{endpoint}", data=data)
        req.raise_for_status()

        return req.json()

    def get_all_studies_metadata(self):
        """Get all studies"""

        return self._getBG('metadata/study')        

    def get_all_tables_metadata(self): 
        """Get all tables"""

        return self._getBG('metadata/table')

    def get_table_metadata(self, table): 
        """Get metadata for a single table"""

        return self._getBG(f'metadata/table/{table}')

    def get_column_metadata(self, table, column): 
        """Get metadata for a single column"""

        return self._getBG(f'metadata/table/{table}/column/{column}')

    def get_all_tissues(self): 
        """Get all available tissues"""

        return self._getBG('metadata/tissue')

    def get_tissue_metadata(self, tissue): 
        return self._getBG(f'metadata/tissue/{tissue}')

    def query(
        self, 
        ids1, 
        ids2=None, 
        columns=None,
        average_columns=False, 
        limit=10000, 
        restriction_join='intersect',
        restriction_bool=None, 
        restriction_gt=None, 
        restriction_lt=None, 
        table='BigGIM_70_v1',
    ): 
        payload = {
            'ids1': ids1, 
            'ids2': ids2, 
            'columns': columns, 
            'average_columns': average_columns, 
            'limit': limit, 
            'restriction_join': restriction_join, 
            'restriction_bool': restriction_bool, 
            'restriction_gt': restriction_gt, 
            'restriction_lt': restriction_lt, 
            'table': table, 
        }

        payload = {
            "ids1":"5111,6996,57697,6815,889,7112,2176,1019,5888,5706",
            "table": table, 
        }

        print(payload)

        out = self._getBG('biggim/query', data=payload)

        return out

