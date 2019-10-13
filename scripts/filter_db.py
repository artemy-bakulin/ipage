import numpy as np


def calculate_matrix(db_profiles, child_unique_genes):
    length = len(db_profiles)
    sums = db_profiles.sum(1)
    matrix = np.memmap('memmapped2.ipage.dat', dtype=np.int16,
                       mode='w+', shape=(length, length))
    for i in range(length):
        matrix[i] = (i + 1) * ((db_profiles[i] > db_profiles).sum(1) / sums[i] < child_unique_genes)
    np.fill_diagonal(matrix, 0)
    matrix -= 1
    matrix = np.transpose(matrix)
    return matrix


def non_redundancy_sort_pre(db_names, db_annotations, db_profiles, child_unique_genes, parent_unique_genes):
    db_sums = np.sum(db_profiles, axis=1)
    db_names = [db_names[i] for i in range(len(db_names)) if db_sums[i] > 6]
    db_annotations = [db_annotations[i] for i in range(len(db_annotations)) if db_sums[i] > 6]
    db_profiles = db_profiles[db_sums > 6]
    db_profiles_memmap = np.memmap('memmapped1.ipage.dat', dtype=np.int16,
                                   mode='w+', shape=np.shape(db_profiles))
    db_profiles_memmap[:] = db_profiles[:]
    matrix = calculate_matrix(db_profiles_memmap, child_unique_genes)

    li = []
    for i in range(len(matrix)):
        indices = [matrix[i] != -1]
        parent_genes_sum = db_profiles[i].sum()
        parent_unique_genes_sum = np.count_nonzero((np.sum(db_profiles[indices], axis=0) - db_profiles[i]) == -1)
        if parent_genes_sum / parent_unique_genes_sum < parent_unique_genes:
            li.append(i)
    db_annotations = [db_annotations[i] for i in range(len(db_annotations)) if i not in li]
    db_profiles = np.array([db_profiles[i] for i in range(len(db_profiles)) if i not in li])
    db_names = [db_names[i] for i in range(len(db_names)) if i not in li]
    return db_names, db_annotations, db_profiles
