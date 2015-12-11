function plotS4grating(grating,lattice)
            strat = grating.stratum{1,1};
	periodX = grating.d21
	periodY = grating.d32
    
    
    figure;
    hold on;
    
    if strcmp(lattice,'square')
        %Define block dimensions and positions:
        %           vertOffset = -periodX;  %Need to shift the block positions up and to the right as you move up in stripe#
        %           horiOffset = -periodY;
        lastX = 0; %marks the position of the end of the last stripe
        for s = 1:length(strat.stripe) %iterate by stripe
            stripe = strat.stripe{1,s};
            stripeEnd = stripe.c1; %top edge of stripe
            lastX = 0; %marks the position of the end of the last block
            for b = 1:length( stripe.block ) %iterate by block
                blockEnd = stripe.block{1,b}.c2; %right edge of current block
                if stripe.block{1,b}.pmt_index == 2 %if it's marked as SU8
                    %                       centerX = ((lastX + blockEnd)/2 + lastY/2) * periodX + horiOffset; %scale by periodX because GDC does it relative to periodicity (period=1) %also shift to right as you go up in y
                    centerX = (lastX + blockEnd)/2 * periodY;
                    widthX = (blockEnd - lastX) * periodY;
                    %                       centerY = (stripeEnd + lastY)/2 * periodY + vertOffset;
                    centerY = (stripeEnd + lastX)/2 * periodX;
                    widthY = (stripeEnd - lastX) * periodX;
                    %setLayerScript = [setLayerScript, sprintf('S:SetLayerPatternRectangle(''Grating'', ''SU8'', {%1.7f,%1.7f}, 0, {%1.7f,%1.7f}) \r\n',centerX,centerY,widthX/2,widthY/2)];
                    
                    rectangle('Position',[centerX-widthX/2, centerY-widthY/2, widthX,widthY]);
                    
                end
                lastX = blockEnd;
            end
            lastX = stripeEnd;
        end
        
        
        %STRIPES PARALLEL TO Y axis
    elseif strcmp(lattice,'hexagonal')
        %Define block dimensions and positions:
        horiOffset = 0;
        vertOffset = 0;
        %horiOffset = -periodX/2;  %Need to shift the block positions up and to the right as you move up in stripe#
        %vertOffset = -periodY/2;
        lastX = 0; %marks the position of the end of the last stripe
        for s = 1:length(strat.stripe) %iterate by stripe
            stripe = strat.stripe{1,s};
            stripeEnd = stripe.c1; %top edge of stripe
            lastY = 0; %marks the position of the end of the last block
            for b = 1:length( stripe.block ) %iterate by block
                blockEnd = stripe.block{1,b}.c2; %right edge of current block
                if stripe.block{1,b}.pmt_index == 2 %if it's marked as SU8
                    %centerY = ((lastY + blockEnd)/2 + lastY/2) * periodX + horiOffset; %scale by periodX because GDC does it relative to periodicity (period=1) %also shift to right as you go up in y
                    centerY = (lastY + blockEnd)/2 * periodY + vertOffset; %scale by periodX because GDC does it relative to periodicity (period=1) %also shift to right as you go up in y
                    %centerX = (lastX + blockEnd)/2 * periodX;
                    widthY = (blockEnd - lastY) * periodY;
                    centerX = (stripeEnd + lastX)/2 * periodX + horiOffset;
                    %centerY = (stripeEnd + lastY)/2 * periodY;
                    widthX = (stripeEnd - lastX) * periodX;
                    %setLayerScript = [setLayerScript, sprintf('S:SetLayerPatternRectangle(''Grating'', ''SU8'', {%1.7f,%1.7f}, 0, {%1.7f,%1.7f}) \r\n',centerX,centerY,widthX/2,widthY/2)];

                    rectangle('Position',[centerX-widthX/2, centerY-widthY/2, widthX,widthY]);
                
                end
                lastX = blockEnd;
            end
            lastX = stripeEnd;
        end
        
    end

end