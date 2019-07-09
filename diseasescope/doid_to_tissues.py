"""Alternate paths to tissues"""

import requests

def get_tissue_from_pmc_w2v(disease, n=5): 
    """Calls GetTissue API (Experimental)
    
    Parameters
    ----------
    disease : str
        Disease name as a string in lower case
    n : int
        Number of tissues to return (in sorted order of most relevant to least)
    """

    out = requests.post(
        "http://secret.ndexbio.org:8085/", 
        data={
            'disease': disease,
            'n': n
        }
    )
    
    if out.status_code == 200: 
        out = out.text.strip()[1:-1]
        out = out.split(',')
        
    else: 
        raise RuntimeError("API call failed! status code: %s" % out.status_code)
    
    return out