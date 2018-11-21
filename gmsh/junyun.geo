 xmin = -0.5; 
 xmax = 0.5;
 ymin = 0; 
 ymax = 1; 
 NX = 16; 
 NY = 16; 
 Point (1) = {xmin, ymin, 0};//设置点
 Point (2) = {xmax, ymin, 0}; 
 Point (3) = {xmax, ymax, 0}; 
 Point (4) = {xmin, ymax, 0}; 
 Line (1) = {1,2};//连线 
 Line (2) = {2,3}; 
 Line (3) = {3,4}; 
 Line (4) = {4,1}; 
 Transfinite Line {1,-3} = NX+1;//在线上布置点 
 Transfinite Line {2,-4} = NY+1; 
 Line Loop (1) = {1,2,3,4};//将封闭的线连成面 
 Plane Surface (1) = {1};//生成面 
 Transfinite Surface {1};//网格生成 
 Recombine Surface {1};//得到结构 
