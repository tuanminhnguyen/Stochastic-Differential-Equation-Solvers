%% ran2 test
num = ran2(-2,1e7);
histogram(num,0:0.01:1)

%% two-point test
t = twopoint(-2,1e5,0.3,0.7,0.2);
histogram(t)

%% exponential test
x = exponential(-2,1e5,2);
histogram(x)

%% standard gaussian tests
x = gauss_box(-2,1e5);
histogram(x)

%% gaussian test
x = gauss_box(-2,1e5,2,3);
histogram(x)

