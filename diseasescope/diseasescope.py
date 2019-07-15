# -*- coding: utf-8 -*-

"""Main module."""

import ndex2
import networkx as nx

from .biggim import BigGIM
from .disgenet import DisGeNet
from .doid_to_genes import doid_to_genes_direct
from .doid_to_tissues import get_tissue_from_pmc_w2v
from .ddot_client import DDOT_Client
from .network import NxNetwork 

class DiseaseScope(object):
    """Implements the DiseaseScope pipline"""

    def __init__(self, doid, disease):
        # self._check_doid(doid)
        self.doid = doid
        self.disease = disease # TODO: Look this up so the disease is not needed


    def get_disease_genes(self, method='biothings'): 
        if method == 'biothings':
            self.genes = doid_to_genes_direct(self.doid)

        elif method == 'disgenet':
            self.disgenet = DisGeNet()
            self.disgenet.query_disease_genes(self.doid, namespace='do')
            self.genes = self.disgenet.get_top_genes()
        
        else: 
            raise ValueError("Invalid method!")

    
    def get_disease_tissues(self, method='pubmed'): 
        if method == 'pubmed': 
            self.tissues = get_tissue_from_pmc_w2v(self.disease)

        elif method == 'biothings': 
            raise NotImplementedError
        
        else: 
            raise ValueError("Invalid method!")


    def expand_gene_set(self, n=250, method='biggim', add_subnetwork=False,**kwargs):
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
            attr = {'heat': {
                i:0 if i in self.genes else 1 for i in self.network.node_names}
            }

            expanded_index = (
                self.network
                    .set_node_attributes(attr)
                    .random_walk('heat', kwargs['alpha'], add_heat=True, heat_name='random walk score')
                    .node_table['random walk score']
                    .sort_values(ascending=False)
                    .iloc[:n]
                    .index
                    .tolist()
            )

            self.expanded_genes = [self.network.node_2_name[i] for i in expanded_index]
        
        elif method == 'heat diffusion': 
            attr = {'heat': {
                i:0 if i in self.genes else 1 for i in self.network.node_names}
            }

            expanded_index = (
                self.network
                    .set_node_attributes(attr)
                    .heat_diffusion('heat', add_heat=True, heat_name='heat diffusion score')
                    .node_table['heat diffusion score']
                    .sort_values(ascending=False)
                    .iloc[:n]
                    .index
                    .tolist()
            )

            self.expanded_genes = [self.network.node_2_name[i] for i in expanded_index]

        else: 
            raise ValueError("Invalid method!")

        if add_subnetwork: 
            self.subnetwork = self.network.subgraph(node_names=self.expanded_genes)

        return self

    def get_network(self, method='biggim',  **kwargs): 
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

            self.edge_table = self.biggim.result_dataframe[['Gene1', 'Gene2', 'mean']] #TODO: Move these column names to attribute
            if kwargs.get('to_network'):
                dG = nx.from_pandas_edgelist(self.edge_table, 'Gene1', 'Gene2', edge_attr=True)
                self.network = NxNetwork(dG)

        elif method == 'ndex':
            if 'uuid' not in kwargs: 
                raise ValueError("uuid is missing from keyword arguments!")

            self.network = NxNetwork().from_ndex(uuid=kwargs['uuid'])

        elif method == 'networkx': 
            self.network = NxNetwork(kwargs['network'])
        
        elif method == 'igraph': 
            from .network import IgNetwork

            self.network = IgNetwork(kwargs['network'])

        else: 
            raise ValueError("Invalid method!")


    def infer_hierarchical_model(
        self, 
        method='clixo-api', 
        method_kwargs=None,
        edge_attr='has_edge',
        ddot_api_location=None,
        get_ontology_object=False,
    ):
        if method_kwargs is None: 
            method_kwargs = {}

        if method == 'clixo-api': 
            cols = self.edge_table.columns
            table = self.edge_table[[cols[0], cols[1], edge_attr]]
            self.ddot_client = DDOT_Client.from_dataframe(table, ddot_api_location=ddot_api_location)
            (self.ddot_client
                .call(**method_kwargs)
                .wait_for_hiview_url()
            )

            self.hiview_url = self.ddot_client.hiview_url

            if get_ontology_object: 
                self.ddot_client.get_output_file(create_object=True)
                self.ontology = self.ddot_client.ontology

        elif method == 'clixo':
            from ddot import Ontology

            ont = Ontology.run_clixo(
                self.edge_table, 
                method_kwargs['alpha'], 
                method_kwargs['beta'], 
                clixo_cmd=method_kwargs['clixo_path'],
            )

            if len(self.edge_table.columns) != 3:
                #TODO
                print("future warning about edge table columns")

            ndex_url, _ = ont.to_ndex(
                method_kwargs['ndex_username'], 
                method_kwargs['ndex_password'], 
                method_kwargs.get('ndex_server', 'http://test.ndexbio.org'), 
                network=self.edge_table, 
                main_feature=self.edge_table.columns[-1]
            )

            uuid = ndex_url.split('/')[-1]
            self.hiview_url = f"http://hiview-test.ucsd.edu/{uuid}?type=test&server=http://dev2.ndexbio.org" 

        else: 
            raise ValueError("Invalid method!")


    def _check_doid(self, doid): 
        raise NotImplementedError

