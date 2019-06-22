import json
import requests
import pandas
import time
from multiprocessing import Pool
import logging

base_url = 'http://biggim.ncats.io/api'

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logging.debug("test")

# =============================================================================
# N = 2 #maximum number of elements to query for
# print(N)
# =============================================================================

def doid_to_genes_and_tissues(doid, direct=True, N=10):
    if direct==True:
        genes = doid_to_genes_direct(doid,N=N)
    else:
        genes = doid_to_genes(doid, N=N)

    tissues = doid_to_tissues(doid, N=N)

    return genes,tissues

def doid_to_genes_direct(doid, N=10): 
    """
    Parameters
    ----------
    doid : list of str
    """

    if not isinstance(doid, list): 
        doid = [doid]

    logging.info("Getting HP ids from DOID")

    #http://disease-ontology.org/
    disease_doid = doid
    genes = doid_to_ncbigene(disease_doid, limit=N)
    if len(genes) == 0:
        raise Exception("No genes related to %s" % (str(disease_doid)))
    else:
        print("Returned %i genes" % (len(genes), ))

    return genes


def doid_to_genes(doid, N=10): 
    """
    Parameters
    ----------
    doid : list of str
    """

    if not isinstance(doid, list): 
        doid = [doid]

    logging.info("Getting HP ids from DOID")

    #http://disease-ontology.org/
    disease_doid = doid
    hpids = doid_to_hp(disease_doid, limit=N)
    if len(hpids) == 0:
        raise Exception("No phenotypes related to %s" % (str(disease_doid)))
    else:
        print("Returned %i phenotypes" % (len(hpids), ))

    logging.info("Geting OMIM from HP ids")
    #https://hpo.jax.org/app/
    omims = hp_to_omim(hpids, limit=N)

    if len(omims) == 0:
        raise Exception("No mendelian diseases(OMIM) related to %s" % (str(hpids)))
    else:
        print("Returned %i mendelian diseases" % (len(omims), ))

    logging.info("Geting genes from OMIM")
    #https://omim.org/
    genes = omim_to_gene(omims, limit=N)
    #genes are ncbi
    if len(genes) == 0:
        raise Exception("No genes associated with %s" % (str(omims),))
    else:
        print("Returned %i genes" % (len(genes,)))

    return genes

def doid_to_tissues(doid, N=10):
    """
    Parameters
    ----------
    doid : list of str
    """

    if not isinstance(doid, list):
        doid = [doid]

    logging.info("Getting HP ids from DOID")
    #http://disease-ontology.org/
    disease_doid = doid
    hpids = doid_to_hp(disease_doid, limit=N)
    if len(hpids) == 0:
        raise Exception("No phenotypes related to %s" % (str(disease_doid)))
    else:
        print("Returned %i phenotypes" % (len(hpids), ))

    logging.info("Getting uberon from HP ids")
    uberontissues = hp_to_uberon(hpids, limit=N)
    if len(uberontissues) == 0:
        raise Exception("No tissues related to %s" % (str(hpids)))
    else:
        print("Returned %i uberon tissues" % (len(uberontissues), ))

    #tissues = uberontissues
    logging.info("Getting brenda tissue names from uberon tissue identifiers")
    tissues = uberon_to_bto(uberontissues, limit=N)
    if len(tissues) == 0:
        raise Exception("No tissues related to %s" % (str(uberontissues)))
    else:
        print("Returned %i brenda tissues" % (len(tissues), ))
         
    return tissues


def query_single_edge(_input, _output, _value):
    """
        Retrieve results using one API call from input to output.
        :param _input: the input prefix, for a list of input prefix in BioThings Explorer, visit: http://biothings.io/explorer/api/v2/metadata/bioentities
        :param _output: the output prefix, for a list of output prefix in BioThings Explorer, visit: http://biothings.io/explorer/api/v2/metadata/bioentities
        :_value: The value of the input
        :return:
    """
    """
    endpoint = 'directinput2output'
    base_url = 'http://biothings.io/explorer/api/v2'
    data = {'output_prefix':_output,
            'input_prefix': _input,
            'input_value': _value,
            'format': 'translator'
    }"""
    #doc = get(endpoint, data=data, base_url=base_url)
    qurl = 'http://biothings.io/explorer/api/v2/directinput2output?'
    qurl += 'input_prefix=%s&output_prefix=%s&input_value=%s&format=translator'
    qurl = qurl % (_input, _output, _value)
    print(qurl)
    try:
        req = requests.get(qurl.strip(), timeout=10)
    except:
        print("timeout")
        return []
    #print("Sent: GET %s?%s" % (req.request.url,req.request.body))
    return req.json()


"""Helpers"""

def el(result):
    try:
        return result['result_list']['edge_list']
    except:
        logging.warning(result)
        return []
    

def _doid_to_ncbigene(doid):

    
    qse = query_single_edge(_input='doid', 
                                _output='ncbigene', 
                                _value=doid
                           )
    results = el(qse)
    ts = set([x['target_id'] for x in results])
    return [x.split(':')[-1] for x in ts if x]

def doid_to_ncbigene(doids, limit=None):
    
    #res = my_pool.map(_doid_to_hp, doids)
    res = []
    
    for doid in doids:
        if limit is None or len(set(res)) < limit:
            res += _doid_to_ncbigene(doid)
    return res

def _doid_to_hp(doid):

    
    qse = query_single_edge(_input='doid', 
                                _output='hp', 
                                _value=doid
                           )
    results = el(qse)
    ts = set([x['target_id'] for x in results])
    return [x.split(':')[-1] for x in ts if x]

def doid_to_hp(doids, limit=None):
    
    #res = my_pool.map(_doid_to_hp, doids)
    res = []
    
    for doid in doids:
        if limit is None or len(set(res)) < limit:
            res += _doid_to_hp(doid)
    return res


def _hp_to_omim(hpid):
    qse = query_single_edge(_input='hp', 
                                _output='omim.disease', 
                                _value=hpid)
    results = el(qse)
    ts = set([x['target_id'] for x in results])
    return [x.split(':')[-1] for x in ts if x]

def hp_to_omim(hpids, limit=None):
    #res = my_pool.map(_hp_to_omim, hpids)
    res = []
    for hpid in hpids:
        if limit is None or len(set(res)) < limit:
            res += _hp_to_omim(hpid)
    return res

def _omim_to_gene(omim):
    qse = query_single_edge(_input='omim.disease', 
                                _output='ncbigene',  
                                _value=omim)
    results = el(qse)
    ts = set([x['target_id'] for x in results])
    return [x.split(':')[-1] for x in ts if x]

def omim_to_gene(omims, limit=None):
    #res = my_pool.map(_omim_to_gene, omims)
    res = []
    for omim in omims:
        if limit is None or len(set(res)) < limit:
            res += _omim_to_gene(omim)
    return res

def _hp_to_uberon(hpid):

    qse = query_single_edge(_input='hp', 
                                _output='uberon', 
                                _value=hpid)
    
    results = el(qse)
    ts = set([x['target_id'] for x in results])
    return [x.split(':')[-1] for x in ts if x]

def hp_to_uberon(hpids, limit=None):
    #res = my_pool.map(_hp_to_uberon,hpids)
    res = []
    
    for hpid in hpids:
        if limit is None or len(set(res)) < limit:
            res += _hp_to_uberon(hpid)
    return res


def uberon_to_bto(ubid, limit=None):
    res = []
    s = ["UBERON:"+y for y in ubid]
    uberon_to_bto = json.load(open('bto_uberon_bg.json'))
    for x in uberon_to_bto:
        if x['uberon_id'] in s and x['bg_label']:
            res += [x['bg_label']]
    return res


# A few helper functions for posting and getting api requests
#a couple of simple helper functions
def post(endpoint, data={}, base_url=base_url):
    req = requests.post('%s/%s' % (base_url,endpoint), data=data)
    req.raise_for_status()
    return req.json()

def get(endpoint, data={}, base_url=base_url):
    req = requests.get('%s/%s' % (base_url,endpoint), data=data)
    try:
        req.raise_for_status()
    except requests.HTTPError as e:
        print("Sent: GET %s?%s" % (req.request.url,req.request.body))
        if e.response.status_code == 400:
            print(e.response)
            jprint(e.response.json())
        raise
    print("Sent: GET %s?%s" % (req.request.url,req.request.body))
    return req.json()

def jprint(dct):
    print(json.dumps(dct, indent=2))



    