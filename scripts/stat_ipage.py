from scipy.stats import hypergeom
import numpy as np
import math
import MI
def hypergeometric(objects_in_bin ,total_size, objects_total, bin_size):
    p_over=math.log(hypergeom.sf(objects_in_bin-1 ,total_size, objects_total, bin_size)+10**(-9),10)
    p_under=math.log(hypergeom.cdf(objects_in_bin,total_size, objects_total, bin_size)+10**(-9),10)
    
    #p_=math.log(hypergeom.pmf(objects_in_bin,total_size, objects_total, bin_size)+10**(-9),10)
    if p_over<p_under:
        p=-p_over
    else:
        p=p_under
    if abs(p)>3:
        return p/(abs(p))*3
    else:
        return p
def get_p_values(go_profile,nbins):
    bin_size=len(go_profile)//nbins
    remain=len(go_profile)-nbins*bin_size
    #remain=0
    p_values=[]
    objects_total=sum(go_profile)
    total_size=len(go_profile)
    objects_in_bin=sum(go_profile[bin_size+remain:bin_size])
    p=hypergeometric(objects_in_bin ,total_size, objects_total, bin_size)
    p_values.append(p)
    for i in range(1,nbins):
        objects_in_bin=sum(go_profile[bin_size*i:bin_size*(i+1)])
        p=hypergeometric(objects_in_bin ,total_size, objects_total, bin_size)
        p_values.append(p)
    return p_values
    
def test_cond_mi(expression_profile,db_profile,abundance_profile, expression_bins, db_bins, abundance_bins,shuffles=1000,max_p=0.005):
    cmi=MI.cond_mut_info(expression_profile,db_profile,abundance_profile, expression_bins, db_bins, abundance_bins)
    max_vectors_over=shuffles*max_p
    expression_profile=expression_profile.copy()
    vectors_over=0
    for i in range(shuffles):
        np.random.shuffle(expression_profile)
        new_cmi=MI.cond_mut_info(expression_profile,db_profile,abundance_profile, expression_bins, db_bins, abundance_bins)
        if cmi<=new_cmi:
            vectors_over+=1
        if vectors_over>=max_vectors_over:
            return False
    return True
    

    