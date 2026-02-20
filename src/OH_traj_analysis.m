clear;clc;
% [x_file,x_path]= uigetfile('*.txt','open_OHdata_file');
% Xname=[x_path,x_file];
% data = importdata(Xname);

% data = load("data_text.mat");
data = load("D:\1Stefan实验室数据\Ca++OH\OH simulation in magnetic trap\trap_20260125\20260125_trajectroy.mat");

m_OH = 17;
% time = data.data2(:,1,:)*1e6;       % mircoseconds(\mu s)
% position = data.data2(:,2:4,:);     % m
% velocity = data.data2(:,5:7,:);     % m/s

time = data.trajectory(:,1,:)*1e6;       % mircoseconds(\mu s)
position = data.trajectory(:,2:4,:);     % m
velocity = data.trajectory(:,5:7,:);     % m/s

L = size(time);
% for i = 1:L(3)
%     plot(time(:,1,i), position(:,1,i))
%     hold on
% end
speed = sqrt(sum(velocity.^2,2));
kinetic_energy = 0.5*m_OH*speed.^2/6.02e26/1.6e-19;   % unit eV
kinetic_energy_muev = kinetic_energy*1e6;                % unit mueV

% inital_time = round(time(1,1));
% final_time = round(time(end,1));

num = 100;

bin = (L(1)-1)/num;
fprintf("The time interval of trajectory is %.2f us\n",bin)
kinetic_energy_min = 0;   % min(kinetic_energy_traget);
kinetic_energy_max = 100; % max(kinetic_energy_traget);
t = kinetic_energy_min:1:kinetic_energy_max;
k=1;
counts = ones(num+1,length(t));
average_kinetic_energy_distribution = ones(num+1,2);
average_kinetic_energy_distribution(:,1) = 0:bin:num*bin;
for i=1:bin:L(1)
    kinetic_energy_traget = kinetic_energy_muev(i,:);
    average_kinetic_energy_distribution(k,2) = mean(kinetic_energy_traget);
    counts(k,:) = histc(kinetic_energy_traget,t);
    k = k + 1;
end

plot(t, counts(end,:))
distribution = [t',counts'];
[a, b] = uiputfile('*.txt', 'Save trajectory TXT file'); % 让用户选择文件名
if a ~= 0   % 检查是否按了“取消”
    filename = fullfile(b, a);  % 拼接完整路径
    save(filename, 'distribution','-ascii');  % 保存变量 trajectory 到 .txt 文件
    fprintf("The file has been save in %s\n",b)
else
    fprintf("No files have been saved!\n")
end

[a, b] = uiputfile('*.txt', 'Save kinetic changing following the time TXT file'); % 让用户选择文件名
if a ~= 0   % 检查是否按了“取消”
    filename = fullfile(b, a);  % 拼接完整路径
    save(filename, 'average_kinetic_energy_distribution','-ascii');  % 保存平均动能随时间的变化到 .txt 文件
    fprintf("The file has been save in %s\n",b)
else
    fprintf("No files have been saved!\n")
end