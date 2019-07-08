# -*- coding: utf-8 -*-

"""Main module."""

from .biggim import BigGIM

class DiseaseScope(object):
    """Implements the DiseaseScope pipline"""

    def __init__(self, doid):
        self._check_doid(doid)
        self.doid = doid

    def get_disease_genes(self, method=''): 
        pass
    
    def get_disease_tissues(self, method=''): 
        pass

    def get_network(self): 
        pass

    def infer_hierarchical_model(self, method='clixo', method_kwargs=None):
        if method_kwargs is None: 
            method_kwargs = {}
        
        pass

    def _check_doid(self, doid): 
        raise NotImplementedError

