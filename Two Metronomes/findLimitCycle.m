function [IC,count] = findLimitCycle(F,timestep)
	time = 0:timestep:1200;
	IC = [0;1]; X = rk42d(F,time,IC);
	count = 1;
	distprev = 0;
	distcurr = 0;
	tolerance = timestep^2;
	len = length(time);
	IC = X(:,end);
	for c = len-1:-1:0
		distprev = distcurr;
		distcurr = ((X(1,c)-IC(1))^2+(X(2,c)-IC(2))^2);
		if(distcurr>=tolerance) 			%when far away from endpoint,
			count = count+1; 				%just increment count of points
		else 
			if(distcurr>distprev) 			%when nearby, check if approaching
				count = count+1; 			%the endpoint or going away from it
			else 							%if approaching, then break 
				break;
			end 
		end
	end
end
