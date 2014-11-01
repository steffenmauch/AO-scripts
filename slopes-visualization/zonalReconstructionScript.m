% simple script which shows zonal reconstruction as well as
% the corresponding slopes
%
% copyright by:
% Steffen Mauch, (c) 10/2014
% email: steffen.mauch (at) gmail.com
%
% You can redistribute it and/or modify it under the terms of the GNU 
% General Public License as published by the 
% Free Software Foundation, version 2.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51
% Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


clc
close all
clear all;

n = 14;             % nb lenses in x/y-direction
f = 4796;           % in units [mu]
pixelSize = 7.4;    % in units [um]
nbPixelsPerLens = 16;
nbPixels = 224;

defaultGrid = 8:nbPixelsPerLens:nbPixels;


xDefault = ones(14,1)*defaultGrid;
yDefault = defaultGrid'*ones(1,14);

dx = zeros(14,14);
dy = zeros(14,14);

dx(7,7) = 0;
dy(7,7) = -1;

%%%%%
[tempX,tempY] = meshgrid(-1:1/n:1-1/n); 
z = exp(-(tempX.^2+tempY.^2) * 1.25 )*10;
intens = tempX;
[dx,dy,intens] = createSpots(z, intens, tempX, tempY, n, n );
%%%%%

Sx = dx*pixelSize/f*nbPixelsPerLens*pixelSize;
Sy = dy*pixelSize/f*nbPixelsPerLens*pixelSize;

W = zonalReconstruction(Sx, Sy, 1);

%W = intgrad2( Sx, Sy );
subplot(2,1,1);
minima = min(min(W));
if( minima < 0 )
   minima = abs(minima);
else
   minima = 0; 
end
surf(W+minima)
colorbar
title('3D plot')
set(gca,'XTickLabel',[0:nbPixelsPerLens*2:nbPixels]*7.4);
set(gca,'YTickLabel',[0:nbPixelsPerLens*2:nbPixels]*7.4);
set(gca,'xtick',[0:2:nbPixelsPerLens]);
set(gca,'ytick',[0:2:nbPixelsPerLens]);
xlabel('x-position in [um]')
ylabel('y-position in [um]')
zlabel('z-height in [um]')

subplot(2,1,2);
quiver(xDefault,yDefault,dx,dy,0)
title('measured slopes')
ylabel('pixel y-axis')
xlabel('pixel x-axis')
grid on
axis([0 224 0 224])
% change grid spacing
set(gca,'xtick',[0:nbPixelsPerLens:nbPixels]);
set(gca,'ytick',[0:nbPixelsPerLens:nbPixels]);

set(gcf, 'PaperPosition', [0 0 20 20]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [20 20]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'test3', 'pdf') %Save figure
