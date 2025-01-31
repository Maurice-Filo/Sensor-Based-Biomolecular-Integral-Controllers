function RevolveFig(Handle_Axis, Handle_Figure, FileName, Angles, Delay, NumColors)
    Handle_Figure.Color = [1, 1, 1];
    [AZ, EL] = view(Handle_Axis);
    for i = 1 : length(Angles)
        view(Handle_Axis, AZ + Angles(i), EL);
        drawnow();
        frame = getframe(Handle_Figure);
        img = frame2im(frame);
        [imgIndexed, colormap] = rgb2ind(img, NumColors); 
        if i == 1
            imwrite(imgIndexed, colormap, FileName, 'gif', 'LoopCount', inf, 'DelayTime', Delay);
        else
            imwrite(imgIndexed, colormap, FileName, 'gif', 'WriteMode', 'append', 'DelayTime', Delay);
        end 
    end
end




