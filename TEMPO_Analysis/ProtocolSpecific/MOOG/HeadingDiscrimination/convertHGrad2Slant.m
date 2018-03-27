function SlantValues = convertHGrad2Slant(varargin)

n = length(varargin);
for i = 1:n
    StimulusGradients(i) = varargin{i};
end


%Convert Horizontal Gradient to Slant Angle
InterocularDistance = 2.8;
ViewingDistance = 30.25;
DegreesVisualAngle = -0.01:0.01;

xLocations = ViewingDistance*tan(DegreesVisualAngle*pi/180);
StimulusDisparities = StimulusGradients*DegreesVisualAngle; %need to add in MeanDisparities here
StimulusDepths = (InterocularDistance./(2*tan(atan(InterocularDistance/(2*ViewingDistance)) - ((StimulusDisparities/2)*pi/180))))  - ViewingDistance;
SlantValues = atan(StimulusDepths/xLocations)*180/pi;

end
