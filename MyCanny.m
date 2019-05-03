function [imx, imy]  = MyCanny(im,sigma,Th,Tl)

im = rgb2gray(im);                                                          %Converting Image into grayscale
[h,w] = size(im);                                                       
% smoothening using a filter

x_gaus = fspecial('gaussian', [1 h],sigma);                                 % X Direction Gaussian
y_gaus = fspecial('gaussian', [w 1],sigma);                                 % Y Direction Gaussian

im_x = imfilter(im,x_gaus);                                                 %Filtering X & Y Direction using seperable filters (Q2)
im_xy = imfilter(im_x,y_gaus);
  
                                                                            % computing gradients
gradX=[-1 0 1;
       -2 0 2; 
       -1 0 1];
          
gradY=[1 2 1;
       0 0 0; 
      -1 -2 -1];
        
GradX= conv2(im_xy,gradX);                                                 %Convoluting both X and Y Gradients
GradY= conv2(im_xy,gradY);
g=abs(GradX)+abs(GradY);                                                   %Obtaining the magnitude of the derivative responses 

                                                                           % computing the angle of gradients
[a,b]=size(GradX);
theta=zeros([a b]);

for i=1:a
      for j=1:b
            if(GradX(i,j)==0)
               theta(i,j)=atan(GradY(i,j)/0.000000000001);
            else
                theta(i,j)=atan(GradY(i,j)/GradX(i,j));
            end
      end
end
 
  theta=theta*(180/pi);
  
  for i=1:a
      for j=1:b
            if(theta(i,j)<0)
                
               theta(i,j)= theta(i,j)-90;
               
                theta(i,j)=abs(theta(i,j));
            end
      end
  end
 
  for i=1:a
      
      for j=1:b
          
          if ((0<theta(i,j))&&(theta(i,j)<22.5))||((157.5<theta(i,j))&&(theta(i,j)<181))
                theta(i,j)=0;
                
          elseif (22.5<theta(i,j))&&(theta(i,j)<67.5)
                 theta(i,j)=45;
                 
          elseif (67.5<theta(i,j))&&(theta(i,j)<112.5)  
                  theta(i,j)=90;
                  
          elseif (112.5<theta(i,j))&&(theta(i,j)<157.5)
                  theta(i,j)=135;
          end
      end
  end 
  
                                                                            %non maximum suppression

st_v=size(im,1);
st_s=size(im,2);

z=zeros(st_v,st_s);
    
for i=2:st_v-1
    
    for j=2:st_s-1 
        
        at=theta(i,j);
        
        if and(-22.5<at,at<=22.5)  
            
            if and(g(i,j)>g(i,j+1),g(i,j)>g(i,j-1)) 
                z(i,j)=g(i,j);
            end
        end
        
        if and(22.5<at,at<=67.5)  
            
            if and(g(i,j)>g(i-1,j-1),g(i,j)>g(i+1,j+1))  
                z(i,j)=g(i,j);
            end
        end
        
        if abs(at)>67.5   
            
            if and(g(i,j)>g(i-1,j),g(i,j)>g(i+1,j))  
                z(i,j)=g(i,j);
            end
            
        end
        
        if and(-67.5<at,at<= -22.5)  
            
            if and(g(i,j)>g(i-1,j+1),g(i,j)>g(i+1,j-1))  
                z(i,j)=g(i,j);
            end
        end
    end
end

g1=im2uint8(g);                                                             %Converting the image into a readble one for our plot 

% grouping edges based on threshold

st_v=size(im,1);
st_s=size(im,2);

hyst_high=g1>=Th;                                                             %hysteris for strong edges

hyst_low=and(g1>=Tl,g1<Th);                                                  %hysteris for weak edges

temp=hyst_high;

place_holder=1;

while place_holder==1  
    place_holder=0;
    for i=2:st_v-1
        
        for j=2:st_s-1
            manjsi=hyst_low(i,j);  
            if manjsi  
                if any(any(hyst_high(i-1:i+1,j-1:j+1)))  
                    hyst_high(i,j)=1;   
                    place_holder=1;   
                    hyst_low(i,j)=0; 
                end
               
            end
        end
    end
end

result=hyst_high;


% edge linking is performed on binary image with edges allready detected.
% for all black pixels, then empty space is connected.

gradient_ok=g1>0.2;  

k_ok=and(theta>-20,theta<20); 

v_ok=and(gradient_ok,k_ok); 

st_v=size(result,1);
st_s=size(result,2);

max_interval=ceil(st_s*0.05); 

for i=1:st_v  
    if result(i,1)==0;  
        connect=1;  
    end
    
    for j=2:st_s
        
        if and(result(i,j-1)==0,result(i,j)==1)  
            link=j-1;  
            
            if (link-connect<max_interval)  
                if all(v_ok(i,connect:link))  
                    result(i,connect:link)=1;    
                end
            end
        end
        
        if and(result(i,j-1)==1,result(i,j)==0)   
            connect=j;  
        end
    end
end

subplot (3, 3, 1),imshow(im);axis image; title('Original Image');
subplot (3, 3, 2),imshow(im_xy);axis image; title('Filtered Image');
subplot (3, 3, 3),imshow(GradX);axis image; title('X-Direction Image');
subplot (3, 3, 4),imshow(GradY);axis image; title('Y- Direction Image');
subplot (3, 3, 5),imshow(g);axis image; title('Magnitude');
subplot (3, 3, 6),imshow(theta);axis image; title('Angle of Pixels');
subplot (3, 3, 7),imshow(g1);axis image; title('Non-Maxima Suppression');
subplot (3,3,8), imshow(result,'InitialMagnification',50);title('Final Edge Detection Image');

end
