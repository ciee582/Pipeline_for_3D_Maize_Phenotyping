function  ShowPoints( points, regions)
%    global Handle_F;
%    if(ishghandle(Handle_F))
%    close(Handle_F);
%    end
   Handle_F = figure('Name','Segmentation Visualization','NumberTitle','off');
   set(gcf,'color','white');

   color = MyGS.MYCOLOR;
   scatter3(points(:,1), points(:,2), points(:,3), 1, [0.0 0.0 0], 'filled');

   hold on;
   for i = 1:length(regions)
       I1 = regions{i};
       
       scatter3(points(I1,1), points(I1,2), points(I1,3), 10, color(i,:), 'filled');
       if i == 1
           scatter3(points(I1,1), points(I1,2), points(I1,3), 10, color(end,:), 'filled');
       end

%        text(points(I1(1),1), points(I1(1),2), points(I1(1),3), num2str(i), 'Color', 'black', 'FontSize', 18);

%        pause(1);
   end
   hold off;

   axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(-180,0);
end

