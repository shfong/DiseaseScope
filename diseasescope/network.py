import pickle
import time
import warnings
from abc import ABC, abstractmethod

import mygene
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import ndex2
import networkx as nx
import scipy
from scipy.sparse import coo_matrix, csc_matrix, csr_matrix, diags, issparse
from scipy.sparse.linalg import expm_multiply


def igraph_adj_matrix(G, weighted=False): 
    source, target, weights = zip(*[(i.source, i.target, i[weighted] if weighted else 1) for i in G.es])

    n_nodes = len(G.vs)
    adj = csr_matrix(coo_matrix((weights, (source, target)), shape=(n_nodes, n_nodes)))

    if not G.is_directed(): 
       adj += adj.T

    return adj

def random_walk_rst(
    F0, 
    A, 
    alpha, 
    normalize=True,  
    axis=1,
    threshold=1e-7,
    max_iter=100, 
    verbose=True
): 
    '''Random walk with restart
    
    Performs random walk with restart on a sparse matrix. If the 
    adjacency matrix is already normalized, this function will 
    work for dense numpy array matrices as well. (set normalize to False)
    
    TODO
    ----
    Update docstring to include normalize change
    Parameters
    ----------
    F0 : scipy.sparse
        Vector or matrix to propagate
    A : scipy.sparse
        Adjacency matrix to propagate over
    alpha : float 
        Restart probability
    threshold : float
        Threshold to consider the propagation has converged
    normalize : bool
        If normalize, the adjacency matrix will be row normalized 
        (divide by the degree)
    axis : int 
        0 or 1. Either row or column normalize
    max_iter: int
        Maximum number of iterations to perform the random walk
    verbose : bool (Deprecated)
        Prints progress (number of iterations and the tolerance)
    '''
    
    counter = 0
    tol = 10
    
    if not issparse(F0) and issparse(A): 
        warnings.warn("Forcing F0 to be sparse") 
        F0 = csr_matrix(F0)
        
 
    if normalize: 
        if issparse(A): 
            A = sparse_normalize(A, axis=axis)
        else: 
            A = dense_normalize(A, axis=axis)

    F_p = F0.copy()
    while tol > threshold: 
        F_t = (1 - alpha)*np.dot(F_p,A) + alpha*F0
        tol = frobenius_norm(F_t - F_p)
        
        F_p = F_t
        counter += 1
        
        if counter > max_iter: 
            warnings.warn('Max iteration reached. Did not converge!')
            
            break
        
    return F_t


def heat_diffusion(laplacian, heat, start=0, end=0.1): 
    """Heat diffusion 
    Iterative matrix multiplication between the graph laplacian and heat
    """

    out_vector=expm_multiply(
        -laplacian, 
        heat, 
        start=start, 
        stop=end, 
        endpoint=True
    )[-1]

    return out_vector


def frobenius_norm(sparse_mat): 
    '''Calculates the frobenius norm of a sparse matrix'''
    
    return np.sqrt(np.power(np.absolute(sparse_mat.data), 2).sum())


def get_common_indices(idx1, idx2):
    '''Gets a set of common index
    
    Take 2 lists and get the intersection of the 
    two lists. Also return the indices needed
    to rearrange each list to get the common index
    '''
    common_idx = np.intersect1d(idx1, idx2)
    
    map1 = dict(zip(list(idx1), range(len(idx1))))
    map2 = dict(zip(list(idx2), range(len(idx2))))
    
    new_idx1 = [map1[i] for i in common_idx]
    new_idx2 = [map2[i] for i in common_idx]
    
    return common_idx, new_idx1, new_idx2


def sparse_normalize(m, axis=0, inplace=False): 
    '''Normalize by one axis
    
    Divide row/column of a sparse matrix by the sum
    of row/column. This implementation does not require the 
    need to create a dense matrix and works directly at
    the coordinates and values of the non-zero elements.
    Parameters
    ----------
    sp_mat : scipy.sparse
        Sparse matrix
    axis : int 
        0/1 (row/column)
        
    Returns
    -------
    mat : scipy.sparse
        row/column normalized sparse matrix
    '''
    if inplace: 
        mat = m
    else:  
        mat = m.copy()
    
    row_index, col_index = mat.nonzero()
    data = mat.data
        
    marginals = np.array(mat.sum(axis=axis)).ravel()

    data = data/marginals[row_index if axis else col_index]
    mat.data = data

    if inplace: 
        return None
    
    return mat


def dense_normalize(m, axis=0, inplace=False): 
    if inplace: 
        mat = m
    else: 
        mat = m.copy()
        
    marginals = np.array(mat.sum(axis=axis))
    marginals[marginals == 0] = 1
    
    mat = mat/marginals
    
    if inplace: 
        return None

    return mat


def calculate_alpha(network, m=-0.02935302, b=0.74842057):
    """Calculate optimal propagation coefficient
    Model from Huang and Carlin et al 2018
    """
    log_edge_count = np.log10(len(network.edges()))
    alpha_val = round(m*log_edge_count+b,3)
    
    if alpha_val <=0:
        # There should never be a case where Alpha >= 1, 
        # as avg node degree will never be negative

        raise ValueError('Alpha <= 0 - Network Edge Count is too high')

    else:
        return alpha_val


class Network(ABC): 
    """Base class for all network classes to inherit from
    
    This base class defines interfacing functionalities that 
    Nbgwas expects. 
    """

    def __init__(self, network=None, node_name = "name"): 
        self.network = network 
        self.node_name = node_name 

        super().__init__() 

    @property
    @abstractmethod
    def adjacency_matrix(self): 
        pass


    @property
    @abstractmethod
    def laplacian_matrix(self): 
        pass


    @abstractmethod
    def add_adjacency_matrix(self): 
        pass 

    @abstractmethod
    def add_laplacian_matrix(self): 
        pass 

    @abstractmethod
    def nodes(self): 
        pass 

    @abstractmethod
    def edges(self): 
        pass

    @abstractmethod
    def subgraph(self):
        pass 

    @abstractmethod
    def get_node_attributes(self): 
        pass

    @abstractmethod 
    def set_node_attributes(self, attr_map, namespace="nodeids"): 
        """set node attributes 
        attr_map is a dictionary of dictionaries
        TODO
        ----
        - Need to make sure the required attributes are created (there are 3 of them)
        """

        pass 

    @property
    @abstractmethod
    def node_ids(self): 
        pass

    @abstractmethod 
    def set_node_names(self, attr=None): 
        pass

    # @property 
    # @abstractmethod 
    # def node_names(self): 
    #     pass

    def convert_node_names(
        self, 
        attribute="name", 
        current="entrezgene", 
        to="symbol", 
        rerun_query=True,
        use_key_for_missing=False, 
        write_to_node_table=True, 
        **kwargs,
    ): 

        """Maps network node names using mygene.info"""

        mg = mygene.MyGeneInfo()
        node_attributes = self.get_node_attributes()

        attr = [v[attribute] for k,v in node_attributes.items()]

        query_result = mg.querymany(
            attr, 
            scopes=current, 
            field=to,
            as_dataframe=True, 
            returnall=True, 
            **kwargs, 
        )

        gene_map = query_result['out'][to].to_dict()

        missing = query_result['missing']
        if missing: 
            if rerun_query: 
                sec_query_df = mg.getgenes(
                    query_result['missing'], 
                    fields='%s,%s' % (current, to),
                    as_dataframe=True
                )

                missing = sec_query_df.loc[sec_query_df['notfound'] == True].index

                gene_map.update(sec_query_df[to].to_dict())

            if len(missing) != 0: 
                warnings.warn('%s nodes cannot be converted. Their original name will be kept!' % len(missing))

                for i in missing: 
                    gene_map[i] = i

        if query_result['dup']: 
            warnings.warn("Gene name conversion contains duplicate mappings!")

        change_to = {}
        for k,v in node_attributes.items(): 
            change_to[k] = gene_map.get(v[attribute], k if use_key_for_missing else None)

        self.set_node_attributes({to:change_to}, namespace="nodeids")

        if write_to_node_table: 
            self.refresh_node_table()

        return self

    def map_attr_data(self, data, store=False): 
        """
        Parameter
        ---------
        data : dict
        """
        
        values = [data.get(node, None) for node in self.node_ids]

        if store: 
            self.set_node_attributes({store: dict(zip(self.node_ids, values))})
        else: 
            return values

    @property
    def node_table(self): 
        if not hasattr(self, "_node_table"): 
            self._node_table = pd.DataFrame.from_dict(dict(self.get_node_attributes())).T
            self._node_table = self._node_table.fillna(0)

        return self._node_table

    @node_table.setter 
    def node_table(self, node_table): 
        #TODO: Add Validation code here

        self._node_table = node_table

    @node_table.deleter
    def node_table(self): 
        if hasattr(self, "_node_table"): 
            del self._node_table

    def refresh_node_table(self): 
        del self.node_table

        self.node_table

        return self

    def refresh_node_attributes(self): 
        self.set_node_attributes(self.node_table.to_dict(), namespace="nodeids")

        return self

    def __getstate__(self): 
        return self.__dict__

    def __setstate__(self, state): 
        self.__dict__.update(state)        

    def to_pickle(self, filename): 
        with open(filename, 'wb') as f: 
            pickle.dump(self, f) 


    @classmethod
    def from_pickle(cls, filename): 
        with open(filename, 'rb') as f: 
            obj = pickle.load(f) 

        return obj

    def random_walk(
        self, 
        node_attr, 
        alpha, 
        add_heat=False, 
        heat_name='diffused heat',
        **kwargs
    ): 

        """Perform random walk"""

        if isinstance(node_attr, str):
            node_attr = [node_attr]

        if isinstance(heat_name, str): 
            heat_name = [heat_name]

        if len(node_attr) != len(heat_name): 
            raise ValueError("node_attr and heat_name needs to have the same number of names!")

        heat = self.node_table.loc[list(self.node_ids), node_attr].values.T
        heat = random_walk_rst(heat, self.adjacency_matrix, alpha, **kwargs)
        heat = np.array(heat.todense())
        
        new_attr = {
            name: {k:v for k,v in zip(self.node_ids, row)} 
                for name, row in zip(heat_name, heat)
        }

        if add_heat: 
            self.set_node_attributes(new_attr)
            self.refresh_node_table()

            return self

        return new_attr
            

class NxNetwork(Network): 
    """Internal object to expose networkx functionalities"""

    def __init__(self, network=None, node_name="name"): 
        super().__init__(network, node_name=node_name)

        if network is not None:
            self.set_node_names(attr=node_name)

        else: 
            self.node_names = None


    @property 
    def node_ids(self): 
        return self.network.nodes()


    def set_node_names(self, attr=None): 
        if attr is None: 
            attr = self.node_name

        self.node_name=attr

        nodes = self.network.node.keys() 

        self.node_names = [
            str(self.network.node[n].get(self.node_name, n)) \
                for n in self.network.nodes()
        ]

        self.node_2_name = dict(zip(nodes, self.node_names))
        self.name_2_node = dict(zip(self.node_names, nodes))

        return self

    
    @property
    def adjacency_matrix(self): 
        if not hasattr(self, "_adjacency_matrix"): 
            self.add_adjacency_matrix()

        return self._adjacency_matrix

    @property
    def laplacian_matrix(self): 
        if not hasattr(self, "_laplacian_matrix"): 
            self.add_laplacian_matrix()

        return self._laplacian_matrix


    def add_adjacency_matrix(self, weights=None): 
        self._adjacency_matrix = nx.adjacency_matrix(
            self.network, weight=weights
        )

        return self


    def add_laplacian_matrix(self, weights=None): 
        self._laplacian_matrix = nx.laplacian_matrix(self.network, weight=weights)

        return self  


    def nodes(self): 
        return self.network.nodes()

    
    def edges(self): 
        return self.network.edges()


    def subgraph(self, node_ids=None, node_names=None): 
        if node_names is not None and node_ids is not None: 
            raise ValueError("Expected either node_names or node_ids. Both given.")

        elif node_names is not None: 
            node_ids = [self.name_2_node[n] for n in node_names]

        return NxNetwork(
            network=self.network.subgraph(node_ids), 
            node_name=self.node_name
        )


    def get_node_attributes(self): 
        return self.network.nodes(data=True) # networkx > 2


    def set_node_attributes(self, attr_map, namespace="nodenames"):
        for attr_name, d in attr_map.items():  
            if namespace == "nodenames": 
                d = {self.name_2_node[k]:v for k, v in d.items() if k in self.name_2_node}

            nx.set_node_attributes(
                self.network, 
                d,
                name=attr_name
            ) # updated to networkx2

        self.refresh_node_table()

        return self


    def from_cx(self, file, node_name="name"):
        """Load CX file as network"""

        del self.__dict__

        network = ndex2.create_nice_cx_from_file(file).to_networkx()

        self.__init__(
            network=network, 
            node_name=node_name
        )
        
        return self


    def from_pickle(self, file, node_name="name"):
        """Read networkx pickle file as network"""

        del self.__dict__
        
        self.__init__(
            network = nx.read_gpickle(file),
            node_name = node_name,
        )

        return self


    def from_ndex(
        self,
        uuid="f93f402c-86d4-11e7-a10d-0ac135e8bacf", #PCNet
        node_name="name",
    ):

        del self.__dict__

        network_niceCx = ndex2.create_nice_cx_from_server(
            server='public.ndexbio.org',
            uuid=uuid
        )

        network = network_niceCx.to_networkx()

        self.__init__(
            network=network, 
            node_name=node_name
        )

        return self


    def to_ndex(
        self,
        name="subgraph",
        server="http://test.ndexbio.org",
        username="scratch2",
        password="scratch2"
    ):

        """Uploads graph to NDEx
        Parameters
        ----------
        name : str
            The key in self.graphs that contains the graph
        server: str
            The NDEx server hyperlink
        username : str
            Username of the NDEx account
        password : str
            Password of the NDEx account
        """

        try:
            g = ndex2.create_nice_cx_from_networkx(self.network)
        except KeyError:
            raise KeyError("%s is not in self.graphs dictionary!" % name)

        uuid = g.upload_to(
            server=server,
            username=username,
            password=password
        )

        return uuid

class IgNetwork(Network): 
    """Internal object to expose igraph functionalities"""

    def __init__(self, network=None, node_name="name"): 
        import igraph as ig

        super().__init__(network, node_name=node_name) 

        if network is not None: 
            self.set_node_names(attr=node_name)

        else: 
            self.node_names = None

    @property
    def node_ids(self): 
        return [v.index for v in self.network.vs]

    def set_node_names(self, attr=None): 
        if attr is None: 
            attr = self.node_name 

        self.node_name = attr

        nodes = [v.index for v in self.network.vs]
        if self.node_name in self.network.vs.attributes(): 
            self.node_names = self.network.vs[self.node_name]
        else: 
            self.node_names = nodes

        self.node_names = [str(i) for i in self.node_names]
        
        self.node_2_name = dict(zip(nodes, self.node_names))
        self.name_2_node = dict(zip(self.node_names, nodes))

        return self

    @property
    def adjacency_matrix(self): 
        if not hasattr(self, "_adjacency_matrix"): 
            self.add_adjacency_matrix()

        return self._adjacency_matrix

    @property
    def laplacian_matrix(self): 
        if not hasattr(self, "_laplacian_matrix"): 
            self.add_laplacian_matrix()

        return self._laplacian_matrix

    def add_adjacency_matrix(self, weights=None): 
        self._adjacency_matrix = igraph_adj_matrix(
            self.network, weighted=weights)

        return self

    def add_laplacian_matrix(self, weights=None): 
        if not hasattr(self, "adjacency_matrix"): 
            self.add_adjacency_matrix(weights=weights)

        D = diags(self.adjacency_matrix.sum(axis=1))
        
        #TODO: Need to test this functionality against networkx
        self._laplacian_matrix = D - self.adjacency_matrix

        return self

    def nodes(self): 
        return self.network.vs

    def edges(self): 
        return self.network.es
    
    def subgraph(self, node_ids=None, node_names=None): 
        if node_names is not None and node_ids is not None: 
            raise ValueError("Expected either node_names or node_ids. Both given.")

        elif node_names is not None: 
            node_ids = [self.node_2_name[n] for n in node_names]

        return self.network.subgraph(node_ids)

    def get_node_attributes(self): 
        attr = {}
        for v in self.network.vs:
            #attr[a] = dict([(i.index, i.attributes()) for i in self.network.vs])
            attr[v.index] = v.attributes()

        return attr

    def set_node_attributes(self, attr_map, namespace="nodenames"): 
        for attr_name, d in attr_map.items():
            attr = [None]*len(self.network.vs)
            for ind, v in enumerate(self.network.vs): 
                if namespace == "nodenames": 
                    if v[self.node_name] in d: 
                        attr[ind] = d[v[self.node_name]]

                elif namespace == "nodeids": 
                    if v.index in d: 
                        attr[ind] = d[v.index] 

                else: 
                    raise ValueError("namespace must be nodenames or nodeids")

            self.network.vs[attr_name] = attr

        return self
