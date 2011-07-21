A1 = load('../GMS_rho_0.01_N_750.txt');
A2 = load('../GMS_rho_0.001_N_750.txt');
A3 = load('../GMS_rho_0.0001_N_750.txt');
A4 = load('../GMS_rho_1e-05_N_750.txt');

B1 = load('../branch1');
B2 = load('../branch2');
B3 = load('../branch3');
B4 = load('../branch4');

x1 = B1(:,1); 
x2 = B2(:,1);
x3 = B3(:,1);
x4 = B4(:,1);

y1 = B1(:,4);
y2 = B2(:,4);
y3 = B3(:,4);
y4 = B4(:,4);

int1 = find(y1==min(y1));
int2 = find(x1==max(x1));
stable = 'linewidth',2;
unstable = '0.8';

figure(1)
plot(B1(1:int1,1),B1(1:int1,4),unstable)
hold on
plot(B1(int1:int2,1),B1(int1:int2,4),stable)
plot(B1(int2:end,1),B1(int2:end,4),unstable)

int3 = find(y2==min(y2));
int4 = find(x2==max(x2));

plot(B2(0:int3,1),B2(1:int3,4),unstable)
plot(B2(int3:int4,1),B2(int3:int4,4),stable)
plot(B2(int4:end,1),B2(int4:end,4),unstable)

int5 = find(y3==min(y3));
int6 = find(x3==max(x3));

plot(B3(int5:end,1),B3(int5:end,4),unstable)
plot(B3(int6:int5,1),B3(int6:int5,4),stable)
plot(B3(1:int6,1),B3(1:int6,4),unstable)

int7 = find(y4==min(y4));
int8 = find(x4==max(x4));

plot(B4(int7:end,1),B4(int7:end,4),unstable)
plot(B4(int8:int7,1),B4(int8:int7,4),stable)
plot(B4(1:int8,1),B4(1:int8,4),unstable)
