clear all
clc
close all

load('DATA_7x7')
load('DATA_8x7')
load('DATA_9x7')
load('VV')

%% plot Ref
hold on 
plot(VV(:,1),VV(:,2))

%% 7x7

DATA=DATA_7x7;
hold on
plot(DATA(:,1),10*log10(DATA(:,2)))


%% 8x7

DATA=DATA_8x7;

% XX=flipud(DATA);
% XX(1)=[];
% DATA=[DATA;XX];

plot(DATA(:,1),10*log10(DATA(:,2)))

%% 9x7



DATA=DATA_9x7;

% XX=flipud(DATA);
% XX(1)=[];
% DATA=[DATA;XX];

plot(DATA(:,1),10*log10(DATA(:,2)))
