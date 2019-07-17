"""Alternate paths to tissues"""

import requests
import json

def get_tissue_from_pmc_w2v(disease, n=5, sep="_"): 
    """Calls GetTissue API (Experimental)
    
    Parameters
    ----------
    disease : str
        Disease name as a string in lower case
    n : int
        Number of tissues to return (in sorted order of most relevant to least)
    """

    response = requests.post(
        "http://diseasescope.ucsd.edu:8085/w2v_tissues",
        data={
            'disease': disease,
            'n': n,
            "sep":sep,
        }
    )
    
    if response.status_code == 200: 
        response_dict = json.loads(response.content)['result']

    else: 
        raise RuntimeError("API call failed! status code: %s" % response.status_code)
    

    return response_dict