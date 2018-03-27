function  m_count =  Docount(data, lbound,  rbound, step)

   for ii = 1 : (rbound-lbound)/step
       
       leftbound = lbound + (ii-1)*step;
       rightbound = lbound + ii*step;
       
       for jj = 1: length(data(1,:))
           
           countindex = find(data(:, jj)< rightbound & data(:,jj) >= leftbound );
           countemp(ii, jj) = length(countindex);
       end
       
       temp(ii,1) = (leftbound + rightbound)/2;
       
   end
   
   
   
   m_count = [temp   countemp];