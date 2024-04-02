%% Exercise 1
load("compEx1.mat")
figure(1)
out = pflat(x2D);
plot(out(1,:), out(2,:), '.')
title('Plot of x2D Matrix')
xlabel('X-axis');
ylabel('Y-axis');

figure(2)
out = pflat(x3D);
plot3(out(1,:), out(2,:), out(3,:), '.')
title('Plot of x3D Matrix')
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis')

%% Exercise 2
figure(3)
load("CompEx2.mat")
img = imread('compEx2.JPG');
imagesc(img);
colormap gray
hold on
plot(p1(1, :), p1(2,:), 'g*', 'MarkerSize', 20); 
plot(p2(1, :), p2(2,:), 'g*', 'MarkerSize', 20);
plot(p3(1, :), p3(2,:), 'g*', 'MarkerSize', 20);

% Calculate lines as the cross product of interconnected points
l1 = cross(p1(:,1), p1(:,2));
l2 = cross(p2(:,1), p2(:,2));
l3 = cross(p3(:,1), p3(:,2));

rital(l1)
rital(l2)
rital(l3)




% Intersection
p_intersect = cross(l2, l3);
p_cartesion = pflat(p_intersect);
plot(p_cartesion(1,1), p_cartesion(2,1), 'r+', 'MarkerSize', 20, 'LineWidth', 5)
hold off


distance = abs(dot(p_cartesion, l1)/(abs(sqrt(l1(1)^2 + l1(2)^2)))) % The absolute value is taken as the relative locations of the line and the points can result in a negative value



%% Exercise 3
load("compEx3.mat")

H_1 = [sqrt(3), -1, 1; 1, sqrt(3), 1; 0, 0, 2]; 
H_2 = [1, -1, 1; 1, 1, 0; 0, 0, 1]; 
H_3 = [1, 1, 0; 0, 2, 0; 0, 0, 1]; 
H_4 = [sqrt(3), -1, 1; 1, sqrt(3), 1; 1/4, 1/2, 2]; 

figure(4)
subplot(2,3,1);
plot([startpoints(1,:); endpoints(1,:)],[startpoints(2,:); endpoints(2,:)], 'b-');
title('Original Points');
axis equal;

% Transformation with H_1
startpoints_H1 = pflat(H_1*[startpoints; ones(1,42)]);
endpoints_H1 = pflat(H_1*[endpoints; ones(1,42)]);
subplot(2,3,2); 
plot([startpoints_H1(1,:); endpoints_H1(1,:)], [startpoints_H1(2,:); endpoints_H1(2,:)], 'b-');
title('Transformation H_1');
axis equal;

% Transformation with H_2
startpoints_H2 = pflat(H_2*[startpoints; ones(1,42)]);
endpoints_H2 = pflat(H_2*[endpoints; ones(1,42)]);
subplot(2,3,3); 
plot([startpoints_H2(1,:); endpoints_H2(1,:)], [startpoints_H2(2,:); endpoints_H2(2,:)], 'b-');
title('Transformation H_2');
axis equal;

% Transformation with H_3
startpoints_H3 = pflat(H_3*[startpoints; ones(1,42)]);
endpoints_H3 = pflat(H_3*[endpoints; ones(1,42)]);
subplot(2,3,4);
plot([startpoints_H3(1,:); endpoints_H3(1,:)], [startpoints_H3(2,:); endpoints_H3(2,:)], 'b-');
title('Transformation H_3');
axis equal;

% Transformation with H_4
startpoints_H4 = pflat(H_4*[startpoints; ones(1,42)]);
endpoints_H4 = pflat(H_4*[endpoints; ones(1,42)]);
subplot(2,3,5); 
plot([startpoints_H4(1,:); endpoints_H4(1,:)], [startpoints_H4(2,:); endpoints_H4(2,:)], 'b-');
title('Transformation H_4');
axis equal;

%% Exercise 4
figure(5)
load("compEx4.mat")
image_1=imread("compEx4im1.JPG");
image_2=imread("compEx4im2.JPG");

camera_center1 = pflat(null([R1 t1])); 
camera_center2 = pflat(null([R2 t2])); 

view_dir1 = R1(end, :)./norm(R1(end, :), 2); % Extracts viewing direction and normalizes the vector
view_dir2 = R2(end, :)./norm(R2(end, :), 2);

U_new = pflat(U);

subplot(1, 2, 1);
plot3(U_new(1,:), U_new(2,:), U_new(3,:), '.'); 
hold on; 
quiver3(camera_center1(1), camera_center1(2), camera_center1(3), view_dir1(1), view_dir1(2), view_dir1(3), 20); 
hold off; 
title('View from Camera 1'); 

subplot(1, 2, 2);
plot3(U_new(1,:), U_new(2,:), U_new(3,:), '.'); 
hold on; 
quiver3(camera_center2(1), camera_center2(2), camera_center2(3), view_dir2(1), view_dir2(2), view_dir2(3), 20); 
hold off; 
title('View from Camera 2'); 






