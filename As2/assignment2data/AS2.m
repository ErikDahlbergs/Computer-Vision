% SCRIPT AS2
%% COMPUTER EXERCISE 1
load compEx1data.mat
% CE1.1
% x_{i} = P_{i}X
figure();
hold on
plot3(X(1,:), X(2,:), X(3,:), 'r.');
plotcams(P)
axis equal

% CE1.2
figure();
img = imread(imfiles{1});
imagesc(img);
colormap gray
hold on
x1 = pflat(P{1}*X);
notnull = isfinite(x{1}(1 ,:));
plot(x1(1,notnull),x1(2,notnull), 'r.')

% CE1.3
T1 = [1 0 0 0; 0 4 0 0; 0 0 1 0; 0.1 0.1 0 1];
T2 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 1/16 1/16 0 1];

X2 = pflat(T1*X);
X3 = pflat(T2*X); 

P1 = P;
P2 = P;

for i=1:9 
    P1{i} = P{i} * T1^-1;
    P2{i} = P{i} * T2^-1;
end 

figure();
plot3(X2(1,:), X2(2,:), X2(3,:), 'r.');
hold on
axis equal
plotcams(P1)

figure();
plot3(X3(1,:), X3(2,:), X3(3,:), 'r.');
hold on
axis equal
plotcams(P2)

% CE1.4
% for T1
figure();
img = imread(imfiles{1});
imagesc(img);
colormap gray
hold on
notnull = isfinite(X2(1 ,:));
a = pflat(P1{1}*X2);
plot(a(1,notnull), a(2,notnull), 'r.');

% for T2
figure();
img = imread(imfiles{2});
imagesc(img);
colormap gray
hold on
notnull = isfinite(X3(1 ,:));
a = pflat(P2{1}*X3);
plot(a(1,notnull), a(2,notnull), 'r.');


%% COMPUTER EXERCISE 2
load compEx1data.mat
T1 = [1 0 0 0; 0 4 0 0; 0 0 1 0; 0.1 0.1 0 1];
T2 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 1/16 1/16 0 1];

%P1 = (P{1}*T1^(-1));
%P2 = (P{2}*T2^(-1));
P1 = (P{2}*inv(T1));
P2 = (P{2}*inv(T2));

% [K, R] = rq(P{1});
[Kt1, Rt1] = rq(P1(:,1:3));
[Kt2, Rt2] = rq(P2(:,1:3));

Kt1 = Kt1./Kt1(3, 3); 
Kt2 = Kt2./Kt2(3, 3);

% Kt1 and Kt2 are not the same transformations

%% COMPUTER EXERCISE 3
load compEx3data.mat
cube1 = imread('cube1.JPG');
cube2 = imread('cube2.JPG');

x1 = pflat(x{1});
x2 = pflat(x{2});

% CE3.1

%mean and std for cube 1
N1 = normise(x1);
N2 = normise(x2);
% CE3.2
x_n1 = N1*x1;
x_n2 = N2*x2;

figure()
plot(x_n1(1,:),x_n1(2,:), 'r*')

figure()
plot(x_n2(1,:),x_n2(2,:), 'r*')

[dim, numPts] = size(Xmodel);

M1 = zeros(numPts * dim, 4 * dim + numPts);
M2 = M1; 


for i = 1:numPts
    for j = 1:dim
        X_pt = Xmodel(:, i)';
        idx = (i - 1) * dim + j;

        startIndex = 4 * (j - 1) + 1;
        endIndex = startIndex + 3;

        M1(idx, startIndex:endIndex) = [X_pt, 1];
        M2(idx, startIndex:endIndex) = [X_pt, 1];

        midIndex = 4 * dim + i;
        M1(idx, midIndex) = -x_n1(j, i);
        M2(idx, midIndex) = -x_n2(j, i);
    end
end

[U1, S1, V1] = svd(M1);
[U2, S2, V2] = svd(M2);

minimal_eig1 = min(S1(S1>0));
minimal_eig2 = min(S2(S2>0));

min1 = V1(:, end);
min2 = V2(:, end);

V1_n = norm(M1*min1); 
V2_n = norm(M2*min2); 

P1 =reshape(min1(1:12),[4 3])'; %the camera matrix
P2 =reshape(min2(1:12),[4 3])';
P1 = inv(N1)*P1; %De-normalized
P2 = inv(N2)*P2; %De-normalized
X = [Xmodel;ones(1,size(Xmodel,2))];
x1_proj = pflat(P1*X);
x2_proj = pflat(P2*X);

c{1} = P1;
c{2} = P2;

figure();
imagesc(cube1);
hold on
plot(x1(1,:), x1(2,:), 'b.', markersize = 20)
plot(x1_proj(1,:), x1_proj(2,:), 'r*')

figure();
imagesc(cube2);
hold on
plot(x2(1,:), x2(2,:), 'b.', markersize = 20)
plot(x2_proj(1,:), x2_proj(2,:), 'r*')

figure();
plot3([Xmodel(1,startind );  Xmodel(1,endind )],...
      [Xmodel(2,startind );  Xmodel(2,endind )],...
      [Xmodel(3,startind );  Xmodel(3,endind)],'b-');
hold on
axis equal
plot3(Xmodel(1,:),Xmodel(2,:),Xmodel(3,:),'x')
plotcams(c)

[R_1,Q_1] = rq(P1);
[R_2,Q_2] = rq(P2);

K_1 = R_1./R_1(3,3)
K_2 = R_2./R_2(3,3)

%% COMPUTER EXERCISE 4
im1 = imread("cube1.JPG");
im2 = imread("cube2.JPG");

[f1 d1] = vl_sift(single(rgb2gray(im1)), 'PeakThresh', 1);
[f2 d2] = vl_sift(single(rgb2gray(im2)), 'PeakThresh', 1);

figure();
imagesc(im1);
hold on
vl_plotframe(f1);

figure();
imagesc(im2);
hold on
vl_plotframe(f2);

[matches ,scores] = vl_ubcmatch(d1,d2);

x1 = [f1(1,matches(1,:));f1(2,matches(1,:))]; 
x2 = [f2(1,matches(2,:));f2(2,matches(2,:))];

perm = randperm(size(matches ,2)); figure;
imagesc([im1 im2]);
hold on;
plot([x1(1,perm(1:10)); x2(1,perm(1:10))+size(im1,2)], [x1(2,perm(1:10)); x2(2,perm(1:10))],'-' );
hold off;

% COMPUTER EXERCISE 5
N = size(x1,2);

% Initialize allPoints
allPoints = zeros(N, 4);

% normalized cameras
P1_n = inv(K_1) * P1;
P2_n = inv(K_2) * P2;

x1n = pflat(inv(K_1) * [x1; ones(1,N)]);
x2n = pflat(inv(K_2) * [x2; ones(1,N)]);
Xmodel = [Xmodel; ones(1, size(Xmodel, 2))];

% for i = 1:N
%     A = [x1(1,i) * P1(3,:) - P1(1,:);
%          x1(2,i) * P1(3,:) - P1(2,:);
%          x2(1,i) * P2(3,:) - P2(1,:);
%          x2(2,i) * P2(3,:) - P2(2,:)];
% 
%     % SVD
%     [something, somethingElse, V] = svd(A);
%     X = V(:, end);
%     matchedPoints(i, :) = X' / X(4);
% end

for i = 1:N
    A = [x1n(1,i) * P1_n(3,:) - P1_n(1,:);
         x1n(2,i) * P1_n(3,:) - P1_n(2,:);
         x2n(1,i) * P2_n(3,:) - P2_n(1,:);
         x2n(2,i) * P2_n(3,:) - P2_n(2,:)];

    % SVD
    [something, somethingElse, V] = svd(A);
    X = V(:, end);
    matchedPoints(i, :) = X' / X(4);
end

x1_proj = pflat(P1*matchedPoints');
x2_proj = pflat(P2*matchedPoints');

figure();
subplot(1, 2, 1);
imagesc(im1);
hold on
plot(x1(1,:), x1(2,:), 'b.', markersize=20);
plot(x1_proj(1,:), x1_proj(2,:), 'r*');

subplot(1, 2, 2);
imagesc(im2);
hold on
plot(x2(1,:), x2(2,:), 'b.', markersize=20);
plot(x2_proj(1,:), x2_proj(2,:), 'r*');

% Extract good points (not more than 3 pixels away from target)
good_points = (sqrt(sum((x1(1:2, :)-x1_proj(1:2, :)).^2)) < 3 & ...
               sqrt(sum((x2(1:2, :)-x2_proj(1:2, :)).^2)) < 3);

figure();
subplot(1, 2, 1);
imagesc(im1);
hold on
plot(x1(1,good_points), x1(2,good_points), 'b.', markersize=20);
plot(x1_proj(1,good_points), x1_proj(2,good_points), 'r*');

subplot(1, 2, 2);
imagesc(im2);
hold on
plot(x2(1,good_points), x2(2,good_points), 'b.', markersize=20);
plot(x2_proj(1,good_points), x2_proj(2,good_points), 'r*');
%%
X = matchedPoints(good_points,:);

figure();
plot3 ([Xmodel(1, startind); Xmodel(1, endind)],...
[Xmodel(2, startind); Xmodel(2, endind )],...
[Xmodel(3, startind); Xmodel(3, endind )] , 'b-');
hold on
plot3(X(:,1), X(:,2), X(:,3),'r.', 'Markersize', 2)
plotcams({P1, P2})
axis equal

function N = normise(x)
    % Calculate the mean of the points
    m = mean(x, 2);

    % Calculate the standard deviation of the points
    s = std(x, 0, 2);

    % Construct the normalization matrix
    N = [1/s(1), 0, -m(1)/s(1); 
         0, 1/s(2), -m(2)/s(2); 
         0, 0, 1];
end



