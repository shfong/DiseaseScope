# -*- coding: utf-8 -*-

"""Main module."""

import mygene
import numpy as np
import ndex2
import networkx as nx
import pandas as pd
from pathlib import Path
from scipy.stats import hypergeom

from .biggim import BigGIM
from .disgenet import DisGeNet
from .doid_to_genes import doid_to_genes_direct
from .doid_to_tissues import get_tissue_from_pmc_w2v
from .ddot_client import DDOT_Client
from .network import NxNetwork 


class Genes(list): 
    """List object to embed additional attributes

    Parameters
    ---------- 
    inList : list
        List-like object used to instantiate the list
    scope : str
        Scope of the gene names. Used in `convert_scope` to convert
        between scopes
    source : str
        Additional string to keep track of where the genes came from.
        Str is not explicitly used for any methods.
    """

    def __init__(self, inList, scope, source=None):
        super().__init__(inList)
        self.scope = scope
        self.source = source
        self.name_map = {}
        
    def convert_scope(
        self, 
        newscope, 
        inplace=False, 
        species="human"
    ):
        """Converts a list of genes from one scope to another"""

        name_map = self.create_name_map(newscope, 
             return_dataframe=False, species=species)

        newList = [name_map.get(i, i) for i in self]

        if inplace: 
            return self.__init__(newList, newscope)
        
        else: 
            return Genes(newList, newscope)

    def create_name_map(self, newscope, return_dataframe=False, species="human"):
        """Creates a dictionary mapping items in the list to a new name"""

        mg = mygene.MyGeneInfo()
        out = mg.querymany(
            self, 
            scopes=self.scope, 
            fields=newscope, 
            species=species, 
            as_dataframe=True
        )

        if return_dataframe: 
            return out

        self.name_map = out[newscope].to_dict()

        return self.name_map



class DiseaseScope(object):
    """Implements the DiseaseScope pipline"""

    BASE_FILE = Path(__file__).parents[1]
    DOID_FILE = BASE_FILE / "data" / "doid" / "doid_name_mappings.txt"
    GENESET_FILE = BASE_FILE / "data" / "genesets" / "combined_genesets.txt"
    PCNET_UUID = "f93f402c-86d4-11e7-a10d-0ac135e8bacf"

    def __init__(self, doid, convert_doid=False, doid_file=None):
        if doid_file is not None: 
            self.DOID_FILE = doid_file

        # self._check_doid(doid)
        self.doid = doid

        if convert_doid: 
            self.disease = self.convert_doid(doid)

    
    def __repr__(self): 
        attr_names = ["seed genes", "genes", "tissues", "network", "hiview_url"] 
        attrs = ["genes", "expanded_genes", "tissues", "network", "hiview_url"]

        names_found = []
        for name, attr in zip(attr_names, attrs): 
            if hasattr(self, attr): 
                names_found.append(name)

                #TODO: Add statistics about attributes? 

        line = f"DiseaseScope query \"{self.disease} (DOID: {self.doid})\" with attributes {', '.join(names_found)}"

        return line

    
    def load_doid_mapping(self, doid_mapping_file = None):
        if doid_mapping_file is None: 
            doid_mapping_file = self.DOID_FILE

        if not Path(doid_mapping_file).exists(): 
            raise RuntimeError(
                f"{doid_mapping_file} cannot be found! Try providing the file in the init. The file can be found in https://github.com/shfong/DiseaseScope/blob/master/data/doid/doid_name_mappings.txt."
            )

        doid_name_mapping = {}
        with open(doid_mapping_file) as f: 
            for line in f: 
                if line[0] == "#": 
                    continue

                key, value = line.split("%")
                doid_name_mapping[key.strip()] = value.strip()

        self.doid_name_mapping = doid_name_mapping

        return self


    def convert_doid(self, doid, load_mapping_file=False):
        if load_mapping_file or not hasattr(self, "doid_name_mapping"): 
            self.load_doid_mapping()

        key = f"DOID:{doid}"

        return self.doid_name_mapping[key]


    def set_disease_genes(self, genes, scope): 
        self.genes = Genes(genes, scope)

        return self

    def set_expanded_genes(self, expanded_genes, scope): 
        self.expanded_genes = Genes(expanded_genes, scope)

        return self

    def get_disease_genes(self, method='biothings'): 
        if method == 'biothings':
            raise NotImplementedError(
                """Biothings API links have changed. Updated links have not """
                """yet been implemented. Please use `disgenet` instead."""
            )
            self.genes = Genes(
                doid_to_genes_direct(self.doid), 
                "entrezgene",
                source="biothings",
            ) 

        elif method == 'disgenet':
            self.disgenet = DisGeNet()
            self.disgenet.query_disease_genes(self.doid, namespace='do')
            self.genes = Genes(
                self.disgenet.get_top_genes(), 
                "symbol",
                source="disgenet",
            )
        
        else: 
            raise ValueError("Invalid method!")

        return self

    
    def get_disease_tissues(self, n=10, method='pubmed', sep=' '): 
        if method == 'pubmed': 
            self.tissues = get_tissue_from_pmc_w2v(self.disease, n=n, sep=sep)

        elif method == 'biothings': 
            raise NotImplementedError
        
        else: 
            raise ValueError("Invalid method!")

        return self


    def expand_gene_set(
        self, 
        n=250, 
        tissue_by_rank=0, 
        method='biggim', 
        add_subnetwork=False,
        **kwargs
    ):

        """Expand the seed gene set

        Parameters
        ----------
        n : int
            Number of total genes to receive
        tissue_by_rank : int
            Which tissue to select (0-index) from the `tissues` attribute
        method : str
            The method used to expand the gene set. Valid options include 
            `biggim`, `random_walk`, `heat_diffusion`. 
        add_subnetwork : bool
            If True, a new attribute will be created called `subgraph` which is
            a subgraph of the current network, by selecting all the nodes 
            returned from `expand_gene_set`.
        """

        if method == 'biggim': 
            self.biggim = BigGIM()
            self.biggim.query(
                ','.join(self.genes), 
                tissues=self.tissues[tissue_by_rank]['tissue'],
                wait_for_query=True, 
                wait_time=1_000_000, 
                to_dataframe=True, 
                average_columns=True, 
                limit=kwargs.get("limit",100_000),
            )
        
            self.expanded_genes = Genes(self.biggim.get_result_genes(n=n), "entrezgene")

        elif method == 'random walk':
            attr = {'heat': {
                i:1 if i in self.genes else 0 for i in self.network.node_names}
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

            self.expanded_genes = Genes([self.network.node_2_name[i] for i in expanded_index], self.genes.scope)
        
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

            self.expanded_genes = Genes([self.network.node_2_name[i] for i in expanded_index], self.genes.scope)

        else: 
            raise ValueError("Invalid method!")

        if add_subnetwork: 
            self.subnetwork = self.network.subgraph(node_names=self.expanded_genes)

        return self

    def get_network(self, tissue_by_rank=0, method='biggim',  **kwargs): 
        if method == 'biggim': 
            if not hasattr(self, "biggim"): 
                self.biggim = BigGIM()

            self.biggim.query(
                ','.join(self.expanded_genes), 
                ids2=','.join(self.expanded_genes), 
                tissues=self.tissues[tissue_by_rank]['tissue'], 
                wait_for_query=True, 
                wait_time=1_000_000, 
                to_dataframe=True, 
                average_columns=True,
                limit=kwargs.get("limit", 100_000)
            )

            self.edge_table = self.biggim.result_dataframe[['Gene1', 'Gene2', 'mean']] #TODO: Move these column names to attribute
            # self.assign_scope_to_edge_table(
            #     {'Gene1': self.expanded_genes.scope, 'Gene2': self.expanded_genes.scope}
            # )

            if kwargs.get('to_network'):
                dG = nx.from_pandas_dataframe(self.edge_table, 'Gene1', 'Gene2', edge_attr=True)
                self.network = NxNetwork(dG)

        elif method == 'ndex':
            if 'uuid' not in kwargs: 
                raise ValueError("uuid is missing from keyword arguments!")

            self.network = NxNetwork().from_ndex(
                ndex_server=kwargs.get("ndex_server", "http://public.ndexbio.org"), 
                ndex_username=kwargs.get("ndex_username", None),
                ndex_password=kwargs.get("ndex_password", None),
                uuid=kwargs['uuid']
            )

        elif method == 'networkx': 
            self.network = NxNetwork(kwargs['network'])
        
        elif method == 'igraph': 
            raise NotImplementedError

            from .network import IgNetwork

            self.network = IgNetwork(kwargs['network'])

        else: 
            raise ValueError("Invalid method!")

        return self



    def get_edge_table(self, attribute, weight=None): 
        network = getattr(self, attribute, None)
        if network is None or not isinstance(network, NxNetwork): 
            raise ValueError("Attribute is not valid!")

        self.edge_table = network.add_edge_table(weight = weight)

        return self



    # def assign_scope_to_edge_table(self, scopes, update=False): 
    #     """Assign scope dictionary to gene name columns to edge_table

    #     #TODO: Might be deprecated

    #     Paramters
    #     ---------
    #     scopes : dict
    #         edge_table columns to scope map
    #     """

    #     if not all([scope in self.edge_table.columns for scope in scopes]): 
    #         raise ValueError("keys to scopes are not in edge_table columns.")

    #     if update: 
    #         self.edge_table_scopes.update(scopes)
    #     else: 
    #         self.edge_table_scopes = scopes

    #     return self
    
    def convert_edge_table_names(
        self, 
        columns,
        scope, 
        newscope, 
        keep=False,
        species="human"
    ):
        """Convert gene names in edge_table"""

        if not hasattr(self, "edge_table"): 
            raise ValueError("No edge_table found!")

        genes = set(self.edge_table[columns].values.ravel().tolist())
        mg = mygene.MyGeneInfo()
        out = mg.querymany(
            genes, 
            scopes=scope, 
            fields=newscope, 
            species=species, 
            as_dataframe=True
        )
        name_map = out[newscope].to_dict()

        new_columns = self.edge_table[columns].applymap(lambda x:name_map.get(x,x))

        if keep:
            new_column_names = [f"{i}_{newscope}" for i in new_columns.columns]
            new_columns.columns = new_column_names 

            self.edge_table = pd.concat([self.edge_table, new_columns], axis=1)
            self.assign_scope_to_edge_table(
                {name:newscope for name in new_column_names}, 
                update=True
            )

        else: 
            self.edge_table.loc[:, columns] = new_columns
            self.assign_scope_to_edge_table(
                {name:newscope for name in columns}, 
                update=True
            )

        return self


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

            self.ontology = Ontology.run_clixo(
                self.edge_table, 
                method_kwargs['alpha'], 
                method_kwargs['beta'], 
                clixo_cmd=method_kwargs['clixo_path'],
                clixo_version=method_kwargs['clixo_version'],
            )

            if len(self.edge_table.columns) != 3:
                #TODO
                print("future warning about edge table columns")
            
            if method_kwargs.get("upload_to_ndex", True): 
                self.upload_ontology_to_ndex(
                    method_kwargs['ndex_username'], 
                    method_kwargs['ndex_password'], 
                    ndex_server = method_kwargs.get('ndex_server', "http://test.ndexbio.org"), 
                    main_network =self.edge_table, 
                    main_feature = self.edge_table.columns[-1]
                )

                ndex_url, _ = self.ontology.to_ndex(
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

        return self

    
    def upload_ontology_to_ndex(
        self, 
        ndex_username, 
        ndex_password, 
        main_network=None, 
        main_feature=None,
        ndex_server="http://test.ndexbio.org"
    ): 

        ndex_url, _ = self.ontology.to_ndex(
            ndex_username, 
            ndex_password, 
            ndex_server, 
            network=main_network, 
            main_feature=main_feature,
        )
        uuid = ndex_url.split('/')[-1]
        self.hiview_url = f"http://hiview-test.ucsd.edu/{uuid}?type=test&server=http://dev2.ndexbio.org" 


    def name_ontology_terms(
        self, 
        method="jaccard", 
        threshold=0.1, 
        correct_for_background=False, 
        update_ontology=False,
        to_ndex=False,
        **kwargs
    ):

        """Performs basic enrichment for each term

        """

        if not hasattr(self, "ontology"): 
            raise ValueError("No ontology found! Please run `infer_hiearchical_"
                "model first!")

        with open(self.GENESET_FILE) as f: 
            combined_genesets = {}
            
            for line in f: 
                key, value = line.split(';')
                value = value.split(',')
            
                combined_genesets[key] = value    

        ontology_genes = np.array(self.ontology.genes)

        if correct_for_background:     
            corrected_threhsold = threshold/len(combined_genesets)
        else: 
            corrected_threhsold = threshold

        term_names = []
        for term, gene_index in self.ontology.term_2_gene.items(): 
            genes_in_term = ontology_genes[gene_index]

            scores = []
            for name, geneset in combined_genesets.items(): 
                if method == "jaccard": 
                    j = jaccard(genes_in_term, geneset)
                    scores.append((name, j))

                elif method == "hypergeometric":
                    pval = hypergeometric(
                        ontology_genes, 
                        geneset, 
                        genes_in_term,
                        background=1
                    )

                    scores.append((name, pval))

            scores = sorted(scores, key=lambda x:x[1], 
                reverse=True if method == "jaccard" else False)

            best_score = scores[0]

            if method == "jaccard": 
                pass_threshold = best_score[1] > corrected_threhsold
            else: 
                pass_threshold = best_score[1] < corrected_threhsold

            if pass_threshold: 
                term_names.append((term, *best_score))
    
        new_name_map = dict([(i,f"{i} - {j}") for i,j,*_ in term_names])

        if update_ontology:
            self.ontology = self.ontology.rename(terms=new_name_map, inplace=True)

        if to_ndex: 
            self.upload_ontology_to_ndex(
                kwargs['ndex_username'], 
                kwargs['ndex_password'],
                ndex_server=kwargs.get("ndex_server",'http://test.ndexbio.org'),
                main_network=kwargs.get("main_network", None), 
                main_feature=kwargs.get("main_feature", None),
            )

        return term_names

    def _check_doid(self, doid): 
        raise NotImplementedError


def jaccard(setA, setB):
    """"Jaccard Index"""

    setA, setB = set(setA), set(setB)
    
    return len(setA.intersection(setB)) / len(setA.union(setB))

def hypergeometric(total_genes, gold_standard_genes, term, background=1): 
    """Performs hypergeometric test"""

    M = len(total_genes)
    N = len(term) # Sample size
    n = len(set(total_genes).intersection(gold_standard_genes))
    x = len(set(term).intersection(gold_standard_genes))
    
    rv = hypergeom.sf(
        x-1, M, n, N
    )
    
    return min(rv*background, 1)
