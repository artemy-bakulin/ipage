import numpy as np
import os


def calculate_matrix(sparse, db_sums, child_unique_genes):
    length = len(sparse)
    matrix = np.zeros((length, length), dtype=np.int8)

    for i in range(length):
        matrix[i] = [len(sparse[i] - sparse[j]) / db_sums[i] < child_unique_genes for j in range(length)]

    np.fill_diagonal(matrix, 0)
    matrix = np.transpose(matrix)
    return matrix


def non_redundancy_sort_pre(db_names, db_annotations, db_profiles, min_pathway_length, child_unique_genes):
    db_sums = np.sum(db_profiles, axis=1)
    db_names = [db_names[i] for i in range(len(db_names)) if db_sums[i] > min_pathway_length]
    db_annotations = [db_annotations[i] for i in range(len(db_annotations)) if db_sums[i] > min_pathway_length]
    db_profiles = db_profiles[db_sums > min_pathway_length]
    db_sums = db_sums[db_sums > min_pathway_length]

    length = len(db_names)

    sparse = []
    for profile in db_profiles:
        sparse.append(set(np.where(profile == 1)[0]))

    matrix = calculate_matrix(sparse, db_sums, child_unique_genes)

    delete = set()
    for i in range(matrix.shape[0]):
        indices = np.where(matrix[i, :] != 0)[0]
        if len(indices) > 1:
            delete.add(i)
        elif len(indices) == 1 and indices[0] > i:
            delete.add(indices[0])

    db_annotations = [db_annotations[i] for i in range(length) if i not in delete]
    db_profiles = np.array([db_profiles[i] for i in range(length) if i not in delete])
    db_names = [db_names[i] for i in range(length) if i not in delete]

    return db_names, db_annotations, db_profiles
