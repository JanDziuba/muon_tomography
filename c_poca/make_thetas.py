#!/usr/bin/python

import os
from os.path import abspath, dirname
os.chdir(dirname(abspath(__file__)))

test_thetas_folder = "../data/test_thetas/"
if not os.path.exists(test_thetas_folder):
    os.makedirs(test_thetas_folder)

train_thetas_folder = "../data/train_thetas/"
if not os.path.exists(train_thetas_folder):
    os.makedirs(train_thetas_folder)


train_filenames = os.listdir("../data/training_data/")
train_filenames = [filename for filename in train_filenames if filename.endswith(".csv")]
train_filenames.sort(key=lambda f: int(f.split('.')[0]))
for filename in train_filenames:
    os.system("./cmake-build-debug/c_poca ../data/training_data/" + filename + " ../data/train_thetas/out" + filename) 

test_filenames = os.listdir("../data/test_data/")
test_filenames = [filename for filename in test_filenames if filename.endswith(".csv")]
test_filenames.sort(key=lambda f: int(f.split('.')[0]))
for filename in test_filenames:
    os.system("./cmake-build-debug/c_poca ../data/test_data/" + filename + " ../data/test_thetas/out" + filename)

