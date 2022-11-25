%% 2
clear;
close all;
clc;
path = '/Users/malteschade/Desktop/Seismic Ac/AESM1590/datasets/H_sht.mat';

% shot_pos=3, trace_offs=2, rec.pos=4, cmp=5
data = load(path).H_SHT;
n_shots = width(unique(data(3,:)));
% shot distance 18m (n=90)

data_cr = sorthdr(data, 4);
cr_d = abs(data_cr(2,end-1));
id = unique(data_cr(4,:));
cr_10 = sum(data_cr(4,:)==id(10));
% receiver distance 9m (n=278)

data_cmp = sorthdr(data, 5);
uni_cmp = unique(data_cmp(5,:));
cmp_d = abs(uni_cmp(1)-uni_cmp(2));
disp(width(uni_cmp));

for i = 1:8
   temp1 = i;
   temp2 = sum(data_cmp(5,:)==uni_cmp(i));
   temp3 = sum(data_cmp(5,:)==uni_cmp(end+1-i));
   temp4 = sum(data_cmp(5,:)==uni_cmp(width(uni_cmp)/2-5+i));
   disp([temp1,temp2,temp3,temp4]);
end

uni_cmp = unique(data_cmp(5,:));
pos = find(data_cmp(5,:)==uni_cmp(width(uni_cmp)/2));
disp(data_cmp(2,pos));

plothdr(data_cmp);