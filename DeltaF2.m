function DF = DeltaF2(RawFluo,window_baseline,window_smooth)
%RawFluo is Neurons x Time
%window_baseline is the window for the moving minimum
%window_smooth is the window for the "filtering" which is a rolling average

if nargin<2
    error('no windows size specified')
end

%rounds even window sizes down to next lowest odd number
if mod(window_baseline,2)==0
    window_baseline=window_baseline-1;
end
if mod(window_smooth,2)==0
    window_smooth=window_smooth-1;
end

%calculates the size of the input data set
n=size(RawFluo);

%Calculates the number of elements in before and after the central element
%to incorporate in the moving mean.  Round command is just present to deal
%with the potential problem of division leaving a very small decimal, ie.
%2.000000000001.
halfspace_baseline=round((window_baseline-1)/2);
halfspace_smooth=round((window_smooth-1)/2);

DF=zeros(size(RawFluo));
for i=1:n(2)
    start=max(1,i-halfspace_baseline);
    stop=min(i+halfspace_baseline,n(2));
    FX=zeros(size(RawFluo,1),stop-start);
    counter=1;
    for j=start:stop
        start_smooth=max(1,j-halfspace_smooth);
        stop_smooth=min(j+halfspace_smooth,n(2));
        FX(:,counter)=sum(RawFluo(:,start_smooth:stop_smooth),2)/(stop_smooth-start_smooth+1);
        counter=counter+1;
    end
    F0=min(FX,[],2);
    DF(:,i)=(RawFluo(:,i)-F0)./F0;    
end