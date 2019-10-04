close all
clear all

basedir = 'E:\pcp2ChR2 data\rebound';
cd(basedir)
load('overallData.mat')

mice = unique(pretestData.mouse);


