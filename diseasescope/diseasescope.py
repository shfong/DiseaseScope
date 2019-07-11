# -*- coding: utf-8 -*-

"""Main module."""

import ndex2

from .biggim import BigGIM
from .doid_to_genes import doid_to_genes_direct
from .doid_to_tissues import get_tissue_from_pmc_w2v
from .ddot import DDOT_Client

class DiseaseScope(object):
    """Implements the DiseaseScope pipline"""

    def __init__(self, doid, disease):
        # self._check_doid(doid)
        self.doid = doid
        self.disease = disease # TODO: Look this up so the disease is not needed

    def get_disease_genes(self, method='biothings'): 
        if method == 'biothings':
            self.genes = doid_to_genes_direct(self.doid)
        
        else: 
            raise ValueError("Invalid method!")
    
    def get_disease_tissues(self, method='pubmed'): 
        if method == 'pubmed': 
            self.tissues = get_tissue_from_pmc_w2v(self.disease)
        
        else: 
            raise ValueError("Invalid method!")

    def expand_gene_set(self, n=250, method='biggim', **kwargs):
        if method == 'biggim': 
            self.biggim = BigGIM()
            self.biggim.query(
                ','.join(self.genes), 
                tissues=self.tissues[0],
                wait_for_query=True, 
                wait_time=1_000_000, 
                to_dataframe=True, 
                average_columns=True,
            )
        
            self.expanded_genes = self.biggim.get_result_genes(n=n)

        elif method == 'random walk':
            raise NotImplementedError
        
        elif method == 'heat diffusion': 
            raise NotImplementedError

        else: 
            raise ValueError("Invalid method!")

    def get_network(self, method='biggim', **kwargs): 
        if method == 'biggim': 
            if not hasattr(self, "biggim"): 
                self.biggim = BigGIM()

            self.biggim.query(
                ','.join(self.expanded_genes), 
                ids2=','.join(self.expanded_genes), 
                tissues=self.tissues[0], 
                wait_for_query=True, 
                wait_time=1_000_000, 
                to_dataframe=True, 
                average_columns=True,
            )

            self.network = self.biggim.result_dataframe[['Gene1', 'Gene2', 'mean']]

        elif method == 'ndex':
            if 'uuid' not in kwargs: 
                raise ValueError("uuid is missing from keyword arguments!")
            network = ndex2.create_nice_cx_from_server(
                server='public.ndexbio.org', 
                uuid=kwargs['uuid'], 
            )
            self.network = network.to_networkx()

        else: 
            raise ValueError("Invalid method!")

    def infer_hierarchical_model(
        self, 
        method='clixo-api', 
        method_kwargs=None
    ):
        if method_kwargs is None: 
            method_kwargs = {}

        if method == 'clixo-api': 
            self.ddot_client = DDOT_Client.from_dataframe(self.network)
            (self.ddot_client
                .call(**method_kwargs)
                .wait_for_hiview_url()
            )


        elif method == 'clixo':
            try: 
                import ddot
            except ImportError("Unable to import DDOT"):
                raise

        else: 
            raise ValueError("Invalid method!")


    def _check_doid(self, doid): 
        raise NotImplementedError

