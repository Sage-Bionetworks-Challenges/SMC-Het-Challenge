import numpy as np

# creates an overlapping matrix from predicted and truth co-clustering matrices

def create_om(ccm_pred, ccm_truth):
    mutations_pred = np.zeros((ccm_pred.shape[0], 1), dtype=int)
    mutations_truth = np.zeros((ccm_truth.shape[0], 1), dtype=int)
    cluster_count = 1
    for row in range(ccm_pred.shape[0]):
        if mutations_pred[row, 0] == 0:
            for i in range(ccm_pred.shape[1]):
                if ccm_pred[row, i] == 1:
                    mutations_pred[i, 0] = cluster_count
            cluster_count+=1

    cluster_count = 1
    for row in range(ccm_truth.shape[0]):
        if mutations_truth[row, 0] == 0:
            for i in range(ccm_truth.shape[1]):
                if ccm_truth[row, i] == 1:
                    mutations_truth[i, 0] = cluster_count
            cluster_count+=1

    num_truth_clusters = mutations_truth.max()
    num_pred_clusters = mutations_pred.max()

    om = np.zeros((num_truth_clusters, num_pred_clusters), dtype=int)

    for i in range(mutations_pred.shape[0]):
        om[mutations_truth[i,0]-1][mutations_pred[i,0]-1] += 1

    return om

if __name__ == "__main__":
    ccm_true1 = np.copy(np.matrix([[1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 1, 0, 0, 0, 0, 0], [0, 1, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 1, 1, 1, 1], [0, 0, 0, 1, 1, 1, 1, 1], [0, 0, 0, 1, 1, 1, 1, 1], [0, 0, 0, 1, 1, 1, 1, 1], [0, 0, 0, 1, 1, 1, 1, 1]]))
    ccm_pred1 = np.copy(np.matrix([[1, 1, 0, 0, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 1, 0, 0, 0, 0], [0, 0, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 0, 1, 1, 1, 1, 1]]))

    print create_om(ccm_pred1, ccm_true1)
    ccm_true = np.copy(np.matrix([[1, 1, 0, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0, 0], [0, 0, 1, 1, 1, 0, 0], [0, 0, 1, 1, 1, 0, 0], [0, 0, 1, 1, 1, 0, 0], [0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 1, 1]]))
    ccm_pred = np.copy(np.matrix([[1, 1, 0, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0, 0], [0, 0, 1, 1, 0, 0, 0], [0, 0, 1, 1, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1], [0, 0, 0, 0, 1, 1, 1], [0, 0, 0, 0, 1, 1, 1]]))
    print create_om(ccm_pred, ccm_true)




