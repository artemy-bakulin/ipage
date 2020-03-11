from scipy.stats import hypergeom
import numpy as np
import math
import MI


def hypergeometric(objects_in_bin, total_size, objects_total, bin_size):
    p_over = math.log(hypergeom.sf(objects_in_bin-1, total_size, objects_total, bin_size)+10**(-9), 10)
    p_under = math.log(hypergeom.cdf(objects_in_bin, total_size, objects_total, bin_size)+10**(-9), 10)
    if p_over < p_under:
        p = -p_over
    else:
        p = p_under
    if abs(p) > 3:
        return p/(abs(p))*3
    else:
        return p


def get_p_values(profile, nbins):
    bin_size = len(profile) // nbins
    remain = len(profile) - nbins * bin_size
    p_values = []
    objects_total = sum(profile)
    total_size = len(profile)
    objects_in_bin = sum(profile[:bin_size+remain])
    p = hypergeometric(objects_in_bin, total_size, objects_total, bin_size + remain)
    p_values.append(p)
    for i in range(1, nbins):
        objects_in_bin = sum(profile[bin_size * i:bin_size * (i + 1)])
        p = hypergeometric(objects_in_bin, total_size, objects_total, bin_size)
        p_values.append(p)
    return p_values


def test_cond_mi(expression_profile, db_profile, abundance_profile, expression_bins, db_bins, abundance_bins,
                 shuffles=10000, max_p=0.005, function='cmi'):
    if function == 'cmi':
        cmi = MI.cond_mut_info(expression_profile, db_profile, abundance_profile, expression_bins, db_bins, abundance_bins)
    elif function == 'mi':
        cmi = MI.cond_mut_info(expression_profile, db_profile, expression_bins, db_bins)
    max_vectors_over = shuffles*max_p
    expression_shuffled_profile = expression_profile.copy()
    vectors_over = 0
    cmis = []
    np.random.seed()
    for i in range(shuffles):
        np.random.shuffle(expression_shuffled_profile)

        if function == 'cmi':
            new_cmi = MI.cond_mut_info(expression_shuffled_profile, db_profile, abundance_profile, expression_bins,
                                       db_bins, abundance_bins)
        elif function == 'mi':
            new_cmi = MI.cond_mut_info(expression_shuffled_profile, db_profile, expression_bins, db_bins)

        cmis.append(new_cmi)
        if cmi <= new_cmi:
            vectors_over += 1
        if vectors_over >= max_vectors_over:
            return 0, False
    cmis.append(cmi)
    z_score = (np.sum(cmis)/len(cmis)-1)/np.std(cmis)
    return z_score, True
