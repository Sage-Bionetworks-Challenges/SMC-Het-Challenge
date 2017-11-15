# Find the correct ordering of the different scenarios for metric behaviour
import csv
import numpy as np

text_dir = './scoring_metric_data/text_files/' # directory to save tsv's and csv's to

# list of scenarios for testing the metric behaviour
scenarios = ["Truth", "OneCluster", "NClusterOneLineage", "NClusterTwoLineages", "NClusterCorrectLineage",
                "ParentIsNieceWithChildren", "ParentIsSiblingWithChildren", "ParentIsCousin","ParentIsAunt", "ParentIsGrandparent", "ParentIsSibling",
                "BigExtraTop", "BigExtraMid", "BigExtraCurBot", "BigExtraNewBot",
                "SmallExtraTop", "SmallExtraMid", "SmallExtraNewBot", 'SmallExtraCurBot',
                "SplitClusterMidMultiChild", "SplitClusterMidOneChild", "SplitClusterBotSame", "SplitClusterBotDiff",
                "MergeClusterTop&Mid", "MergeClusterMid&BotMultiChild", "MergeClusterMid&BotOneChild","MergeClusterBot"]

def get_corr_order(in_file, out_file=None):
    '''Find the overall best ordering of the scenarios by taking into account the orderings of multiple stakeholders.
    Calculates the overall ordering using multiple different aggregating methods and returns each one, giving you a choice
    of which to use.
    Methods of calculating the overall correct order: Copeland, Borda

    :param in_file: filepath to a csv file containing stakeholders ordering of each scenario
    :param out_file: file to output the data used to calculate the correct order and the overall correct orderings found
    :return:
    '''
    with open(in_file) as csvfile:
        r = csv.DictReader(csvfile, delimiter=',')
        names = np.copy(r.fieldnames)
        names[0] = 'Truth'

        n_scenario = len(names) # Add one for the truth scenario, minus one for the column with the stakeholders' names

        # for each pair of scenarios find the number of stakeholders that rank the
        # first scenario over the second and vice versa
        comparisons = np.array([[0] * n_scenario] * n_scenario)
        data = list()
        for row in r:
            row['Truth'] = n_scenario
            data.append(row)
            for i in range(n_scenario):
                sc1 = names[i]
                for j in range(i+1, n_scenario): # calculate 'winner' for each pair of scenarios for the given ranking
                    sc2 = names[j]
                    if sc1 == 'Truth':
                        comparisons[i][j] += 1
                    elif sc2 == 'Truth':
                        comparisons[j][i] += 1
                    elif float(row[sc1]) > float(row[sc2]):
                        comparisons[i][j] += 1
                    elif float(row[sc1]) < float(row[sc2]):
                        comparisons[j][i] += 1
                    else:
                        comparisons[i][j] += 0.5
                        comparisons[j][i] += 0.5

        # calculate the Copeland ranking
        copeland = np.sum(comparisons, axis=1) - np.sum(comparisons, axis=0)

        # calculate the Borda ranking
        borda = np.array([0] * n_scenario)
        for row in data:
            borda = np.array([float(row[names[i]]) for i in range(n_scenario)]) + borda

        # calculate the standard deviation in the ranking of each scenario
        std_dev = [np.std([float(row[sc]) for row in data]) for sc in names]

        if out_file is not None:
            with open(out_file, 'wb') as f:
                w = csv.writer(f, delimiter=',')
                w.writerow(['Case'] + [row['Ordering'] for row in data] + ['Copeland', 'Borda', 'Std Dev'])
                for i in range(n_scenario):
                    w.writerow([names[i]] + [row[names[i]] for row in data] + [copeland[i]] + [borda[i]] + [std_dev[i]])

        return names, copeland, borda
