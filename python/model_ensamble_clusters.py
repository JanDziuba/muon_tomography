# ---------------------------------
# --- RANDOM FORREST BIN - 100 ----
# --- RANDOM FORREST PRED - 150 ---
# ---------------------------------
# -------- theta_vals = 10 --------
# -- clusters [count x av_theta] --	| 2						| 5						| 10					| 40
# ------ 0.2 < z axis < 0.2 -------
# ----- fill 0 if count < 10 ------
# ------------ count in -----------
# ------------- neigh -------------
# Accuracy bin:						| 0.966592875			| 0.96903975			| 0.96955025			| 0.969550625
# Accuracy pred:					| 0.892982375			| 0.89498375			| 0.895364875			| 0.89544375
# score:
# (offset)	14:						| 0.04200295091397502	| 0.04195688082531596	| 0.041679132835804274	| 0.041280252179147166
#
# input:	output_theta_vals7		| 
# output:							| pred_ensamble_bin11	| pred_ensamble_bin12	| pred_ensamble_bin13	| pred_ensamble_bin14
#									| pred_ensamble13		| pred_ensamble14		| pred_ensamble15		| pred_ensamble16
#
# scores:							| 14_neigh.csv			| 14_1_neigh.csv		| 14_2_neigh.csv		| 14_3_neigh.csv
# ---------------------------------
# ---------------------------------
# --- RANDOM FORREST BIN - 100 ----
# --- RANDOM FORREST PRED - 150 ---
# ---------------------------------
# -------- theta_vals = 10 --------
# -- clusters [count x av_theta] --	| 10					| 2
# ------ 0.2 < z axis < 0.2 -------
# ----- fill 0 if count < 10 ------
# ------------ count in -----------
# ----------- neigh full ----------
# Accuracy bin:						| 0.969025875			| 
# Accuracy pred:					| 0.8960345				| 
# score:
# (offset)	15:						| 0.041436279971244445	| 
#
# input:	output_theta_vals7		| 
# output:							| pred_ensamble_bin15	| pred_ensamble_bin16
#									| pred_ensamble17		| pred_ensamble18
#
# scores:							| 15_neigh.csv			| 15_1_neigh.csv
# ---------------------------------



import pandas as pd
import numpy as np

import re
import os
import glob
import gc


all_files_data = glob.glob(os.path.join("output/output_theta_vals7", "*.csv"))
all_files_data.sort(key=lambda f: int(re.sub('\D', '', f)))
train_data = pd.concat((pd.read_csv(f, header=None) for f in all_files_data[:15000]), ignore_index=True)
test_data = pd.concat((pd.read_csv(f, header=None) for f in all_files_data[15000:]), ignore_index=True)

all_files_target = glob.glob(os.path.join("/home/wnukers/Documents/phd/programs/geometry3/data", "*_ground"))
all_files_target.sort(key=lambda f: int(re.sub('\D', '', f)))
target_train = pd.concat((pd.read_csv(f, header=None, sep=' ', dtype="int8").stack() for f in all_files_target[:15000]), ignore_index=True)
target_test = pd.concat((pd.read_csv(f, header=None, sep=' ', dtype="int8").stack() for f in all_files_target[15000:]), ignore_index=True)


train_d = np.array(train_data[train_data.columns[1:10]]).flatten().reshape(-1)
train_d = train_d.reshape(15000,14400)
train_count = np.array(train_data[train_data.columns[0]]).flatten().reshape(-1)
train_count = train_count.reshape(15000,1600)

test_d = np.array(test_data[test_data.columns[1:10]]).flatten().reshape(-1)
test_d = test_d.reshape(5000,14400)
test_count = np.array(test_data[test_data.columns[0]]).flatten().reshape(-1)
test_count = test_count.reshape(5000,1600)

data = np.concatenate([train_d,test_d])
data = data.reshape(20000,14400)
data_c = np.concatenate([train_count,test_count])
data_c = data_c.reshape(20000,1600)

means = np.array([[np.mean(ar_c),np.mean(ar_d)] for ar_c,ar_d in zip(data_c,data)])


del all_files_data, all_files_target, train_d, test_d, data, data_c, train_count, test_count
gc.collect()


# Cluster means
from sklearn.cluster import KMeans

n_of_clusters = 2
kmeans = KMeans(n_of_clusters)
kmeans.fit(means)
clusters = kmeans.fit_predict(means)

cluster_train = np.repeat(clusters[:15000], 1600)
cluster_test = np.repeat(clusters[15000:], 1600)

train_data_clusters = [train_data[cluster_train==i] for i in range(n_of_clusters)]
train_target_clusters = [target_train[cluster_train==i] for i in range(n_of_clusters)]
test_data_clusters = [test_data[cluster_test==i] for i in range(n_of_clusters)]

del kmeans, cluster_train, clusters
gc.collect()


# Random forrest - binary class
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

dec_tree_bin_cluster_preds = []
for c in range(n_of_clusters):
	bin_target_train = np.where(train_target_clusters[c] == 0, False, True)
	crf = RandomForestClassifier(n_estimators=100, n_jobs=20)
	crf.fit(train_data_clusters[c], bin_target_train)
	dec_tree_bin_cluster_preds += [crf.predict(test_data_clusters[c])]

	print(c)

	del crf, bin_target_train
	gc.collect()


dec_tree_bin_preds = []
idxs = [0 for i in range(n_of_clusters)]
for c in cluster_test:
	dec_tree_bin_preds += [dec_tree_bin_cluster_preds[c][idxs[c]]]
	idxs[c] += 1


off = np.zeros((40,40), dtype='bool')
off[5:35,5:35] = 1
off_test = np.concatenate([off for i in range(5000)], axis=None)

dec_tree_bin_preds = np.array(dec_tree_bin_preds)
dec_tree_bin_preds *= off_test

dec_bin_preds = dec_tree_bin_preds.reshape((5000,40,40)).astype(bool)
# target_pred_neigh = np.zeros((5000,40,40))

neigh_idx = 0
flag = True
while flag:
	print(neigh_idx)

	flag = False
	target_pred_neigh = np.zeros((5000,40,40))
	for idx,grid in enumerate(dec_bin_preds):
		for x in range(6,35):
			for y in range(6,35):
				neigh = np.sum([grid[x-1,y-1], grid[x,y-1], grid[x+1,y-1], grid[x-1,y], grid[x+1,y], grid[x-1,y+1], grid[x,y+1], grid[x+1,y+1]])

				if grid[x,y] == 0:
					if neigh > 5:
						target_pred_neigh[idx][x,y] = 1
						flag = True
				else:
					target_pred_neigh[idx][x,y] = 1

	dec_bin_preds = target_pred_neigh
	neigh_idx += 1

for idx,grid in enumerate(dec_bin_preds):
	for x in range(6,35):
		for y in range(6,35):
			neigh = np.sum([grid[x-1,y-1], grid[x,y-1], grid[x+1,y-1], grid[x-1,y], grid[x+1,y], grid[x-1,y+1], grid[x,y+1], grid[x+1,y+1]])

			if grid[x,y] == 1 and neigh < 3:
				target_pred_neigh[idx][x,y] = 0


target_pred_neigh = target_pred_neigh.flatten().astype(bool)
dec_tree_bin_cluster_preds_neigh = []
for c in range(n_of_clusters):
	dec_tree_bin_cluster_preds_neigh += [target_pred_neigh[cluster_test==c].astype(bool)]


bin_test_target = np.where(target_test == 0, False, True)
print("Accuracy:", accuracy_score(bin_test_target, dec_tree_bin_preds))
print("Accuracy:", accuracy_score(bin_test_target, target_pred_neigh))

idx = 15000
for grid in target_pred_neigh.reshape(5000,40,40):
	gr_name = "pred/pred_ensamble_bin16/" + str(idx) + "_pred"
	np.savetxt(gr_name, grid, fmt='%i')
	idx += 1

del target_pred_neigh, dec_tree_bin_preds, bin_test_target, dec_bin_preds, dec_tree_bin_cluster_preds, off
gc.collect()


# Random forrest - multiclass
dec_tree_class_cluster_preds = []
for c in range(n_of_clusters):
	crf = RandomForestClassifier(n_estimators=150, n_jobs=20)
	crf.fit(train_data_clusters[c][train_target_clusters[c]!=0], train_target_clusters[c][train_target_clusters[c]!=0])
	dec_tree_class_cluster_preds += [crf.predict(test_data_clusters[c][dec_tree_bin_cluster_preds_neigh[c]])]

	print(c)

	del crf
	gc.collect()


idxs_bin = [0 for i in range(n_of_clusters)]
idxs_class = [0 for i in range(n_of_clusters)]
dec_tree_preds = []
for c in cluster_test:
	if dec_tree_bin_cluster_preds_neigh[c][idxs_bin[c]]:
		dec_tree_preds += [dec_tree_class_cluster_preds[c][idxs_class[c]]]
		idxs_class[c] += 1
	else:
		dec_tree_preds += [0]
	idxs_bin[c] += 1


dec_tree_preds = np.array(dec_tree_preds)
dec_tree_preds *= off_test

print("Accuracy:", accuracy_score(target_test, dec_tree_preds))


idx = 15000
for grid in dec_tree_preds.reshape(5000,40,40):
	gr_name = "pred/pred_ensamble18/" + str(idx) + "_pred"
	np.savetxt(gr_name, grid, fmt='%i')
	idx += 1


# Evaluation
def aIOU(predictions, gt):
    classes = np.unique(np.concatenate((np.unique(predictions), np.unique(gt))))
    IoUs = np.zeros(len(classes) - 1)
    i = 0
    for cls in classes:
        if cls != 0:
            preds_tmp = predictions == cls
            gt_tmp = gt == cls
            IoUs[i] = np.sum(np.logical_and(preds_tmp, gt_tmp)) / np.sum(np.logical_or(preds_tmp, gt_tmp))
            i += 1
    return np.mean(IoUs)

dec_tree_eval = dec_tree_preds.reshape((5000,1600))
target_eval = np.array(target_test).reshape((5000,1600))

score = np.mean([aIOU(dec_tree_eval[i], target_eval[i]) for i in range(len(dec_tree_eval))])
print(score)

# Output scores
output = np.array([dec_tree_preds.reshape(5000,40,40)[i].T.reshape(1600) for i in range(5000)])
np.savetxt("scores/15_1_neigh.csv", output, fmt='%i', delimiter=',')
