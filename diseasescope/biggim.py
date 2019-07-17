"""BigGIM Client"""
import requests
import logging
from time import time, sleep
import pandas as pd
import numpy as np
from numpy import inf

logger = logging.getLogger()
logger.setLevel(logging.INFO)

class BigGIM(object):
    """Client for BigGIM
    
    Inspired by this repository: 
    https://github.com/NCATS-Tangerine/BigGIM_APIWrapper/blob/master/app.py
    that wraps the BigGIM API.

    As well as Wrapper code by John Earls 
    TODO: Include links
    """

    URL = "http://biggim.ncats.io/api"
    API_DOC = "http://biggim.ncats.io/api/#!/biggim/get_interactions_query"

    def __init__(self, query_id=None):
        self.query_id = query_id
        
        if query_id is not None: 
            self.check_query_status()

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

    def get_all_columns(self, table):
        """Get all available columns in a table"""

        tmp = self.get_table_metadata(table)
        columns = [j['name'] for j in tmp['columns']]

        return columns


    def get_column_metadata(self, table, column): 
        """Get metadata for a single column"""

        return self._getBG(f'metadata/table/{table}/column/{column}')

    def get_all_tissues(self): 
        """Get all available tissues"""

        return self._getBG('metadata/tissue')

    def get_tissue_metadata(self, tissue): 
        return self._getBG(f'metadata/tissue/{tissue}')

    def query_columns(self, tissues, table='BigGIM_70_v1'):
        """Return the relevant columns from BigGIM to the query tissue

        Parameters
        ----------
        tissues : str or iterable
            Query or queries to retrieve relevant columns from BigGIM. If
            tissues a str, it is assumed to be a single query and will be
            converted to a list with a single item
        table : str
            Name of the BigGIM table to look up
        """

        #TODO: Refactor this

        columns = []
        try:
            # query
            md = self.get_tissue_metadata(tissues)
            for st in md.get('substudies', []): # Allow for nothing found
                for column in st.get('columns',[]):
                    #filter table
                    if column.get('table', {}).get('name', None) == table:
                        #filter r(spearman)
    #                    if column.get('interactions_type',None) == 'Spearman Rank Correlation Coefficient':
                            #get the name of the columns
                            if column.get('name', None) is not None:
                                columns.append(column.get('name'))
        except requests.HTTPError as e:
            #if 404 error, it could not find the tissue, otherwise something bad happend
            if e.args[0].find('404') == -1:
                raise
        
        return columns



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
        tissues=None, 
        wait_for_query=False, 
        wait_time = 1, 
        to_dataframe = False, 
    ):
        """
        Thin wrapper around BigGIM query API 
        
        Parameters
        ----------
        ids1: str
            Comma separated list of gene names
        ids2: str 
            Optional second set of genes to cross with ids1. Otherwise, all 
            genepairs containing ids1 will be considered
        columns : str
            Which BigGIM table column to consider.
        average_columns : Bool
            If True, the data are averaged (ignoring nans). This is done on the
            client side instead of the server (API side)
        limit : int
            Maximum number of genepair to consider
        restriction_join: bool 
            TODO: ???
        restriction_bool 
        restriction_gt 
        restriction_lt 
        table : str
            Name of the BigGIM table to use
        tissues : str
            columns take precedence over tissues. If tissues is given, the 
            columns tha correspond to the tissue will be returned.
        wait_for_query : Bool
            If true, the method call will not return until the query is complete
        wait_time : int
            Maximum number of seconds to wait. After the time limit is reach, 
            the current status of the query will be returned.
        to_dataframe : Bool
            If true, the data from BigGIM will be returned as a dataframe. 
            Otherwise, the links to the data will be located in `request_uri` 
            attribute if the query was completed successfully.
        """

        if columns is not None: 
            validated_columns = columns

        elif tissues is not None:
            validated_columns = self.query_columns(tissues, table=table)
            logging.info(f"Found the following columns in BigGIM: {validated_columns}")
            
        else: 
            validated_columns = None # This is strip columns from the input, resulting in all columns being returned

        validated_columns = ','.join(validated_columns)

        # average_columns is left off to be averaged later 
        kwargs = {
            'ids1': ids1, 
            'ids2': ids2, 
            'columns': validated_columns, 
            'limit': limit, 
            'restriction_join': restriction_join, 
            'restriction_bool': restriction_bool, 
            'restriction_gt': restriction_gt, 
            'restriction_lt': restriction_lt, 
            'table': table, 
        }

        payload = {k:v for k,v in kwargs.items() if k is not None}

        logging.debug(f"Payload: {payload}")
        self.query_payload = payload

        try: 
            self.response = self._getBG('biggim/query', data=payload)
            logging.debug(f"Response received: {self.response}")

            if self.response.get('status', None) != 'submitted' or not self.response.get('request_id', None):
                raise RuntimeError("BigGIM query failed with the following response: %s" % str(self.response))

        except requests.HTTPError as e: 
            logging.info("Request to BigGIM failed") 
            logging.debug(self.response.json())

            raise e

        logging.info("Successfully submitted query")
        logging.debug(f"Response from BigGIM: {self.response}")
        self.query_id = self.response['request_id']

        if wait_for_query: 
            self.wait_for_query(timelimit=wait_time)

        if to_dataframe: 
            return self.convert_result_to_dataframe(average_columns=average_columns)

        return self

    def wait_for_query(self, timelimit=None):
        """Wait for the query to complete

        Paramters
        ---------
        timelimit : int
            Maximum time (in seconds) to wait for the query. If the query is not
            completed yet, a RuntimeError will be raised.
        """
        start_time = time()

        if timelimit is None: 
            timelimit = inf

        if not hasattr(self, "query_id"):
            raise ValueError("No query_id found! Please submit a query first.")

        sleep_time = 1
        while True: 
            self.check_query_status()
            if self.status != 'running': 
                break

            elapsed_time = time() - start_time
            if elapsed_time > timelimit:
                raise RuntimeError("Query over the time limit")
            
            sleep(sleep_time)
            sleep_time += 1

        return self

    def convert_result_to_dataframe(self, average_columns=False): 
        """Retrieve the results and convert it into a single dataframe
    
        The first column is dropped right now (not sure why) TODO    
        """

        if not hasattr(self, "response"): 
            raise ValueError("No response attribute found. Try submiting a query.")
        
        if 'request_uri' not in self.response:
            raise ValueError("Query either did not complete yet or did not complete successfully.")
        
        self.request_uri = self.response['request_uri']
        df = pd.concat(
            [pd.read_csv(i) for i in self.request_uri]
        ).iloc[:, 1:]

        df = df.astype({'Gene1': str, 'Gene2': str})

        if average_columns: 
            df['mean'] = df.iloc[:, 2:].apply(np.nanmean, axis=1)
        
        self.result_dataframe = df

        return self

    def check_query_status(self):
        if not hasattr(self, "query_id"): 
            raise ValueError("No query_id found. Try submitting a query.")

        self.response = self._getBG(f'biggim/status/{self.query_id}')
        self.status = self.response['status']

        if self.status != 'running': 
            logging.info("Query has complete!")

            self.request_uri = self.response.get('request_uri', None)

        return self

    def get_result_genes(self, column='mean', n=None, sort_ascending=False):
        """Get genes ranked by the BigGIM interaction results

        TODO: Make aggregation method visible
        - Right now, the column 'mean' is hard coded

        Paramters
        ---------
        n : int 
            Top n genes, in addition to the query, to be returned. If n is None,
            all genes will be returned. 
        sort_ascending : bool 
            If True, the genes will be sorted in ascending order. 
        """
        if not hasattr(self, "result_dataframe"): 
            raise ValueError("No result_dataframe found. Try submitting "
                "a query and converting to dataframe") 
        
        if n is None: 
            n = inf

        df  = (self.result_dataframe
            .sort_values(by='mean',ascending=sort_ascending)
            .reset_index(drop=True)
        )

        # del df['index'] #drop should make this not necessary

        expanded_genes = set(self.query_payload['ids1'].split(','))
        for _, row in df.iterrows():
            if len(expanded_genes) >= n:
                break

            expanded_genes = expanded_genes.union([int(row['Gene1']), int(row['Gene2'])]) 

        expanded_genes = sorted([str(int(x)) for x in expanded_genes])

        return expanded_genes


        