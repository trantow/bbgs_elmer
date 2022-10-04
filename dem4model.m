% Make bering surface for model
% Need a rectangular box as a DEM -> use .dtm with 0's where dtm is not
% calculated inside the box
% Thomas Trantow
% 07/17/2013
%--------------------------------------------------------------------------
% Input:
%        dtmName: Path to dtm file from kriging output
%        outputname: Name of ouput file ready to use in Elmer
%        Optional - Resolution of DEM (Default 200m)
%
% Note: One will need to change the x and y values if the DEM is other than
%       Bering Glacier or The Bering-Bagley System. Values for the DEM Box
%       for the BBGS and Bering only are given in the code. BBGS is default
%--------------------------------------------------------------------------

function dem4model(dtmName, outputName, varargin)
%DEM resolution (default 200 meters)
if nargin ==2
  res = 200; %meters
elseif nargin ==3
    res = varargin{1};
else
    error('Invalid Number of Arguments: Please give three arguments: (1) The file name of the dtm file for conversion, (2) The name of the output file ready for Elmer and (3) Optional: The DEM resolution (default is 200 meters if nothing is entered)')
        
end


% DEM box
% Bering
%x = (336000:res:444800);
%y = (6660000:res:6719800);
% Bering - Bagley (with glims)
x = (336000:res:508000);
y = (6660000:res:6724000);

% Unikrg calculation of DEM
if isstr(dtmName)
    DTM = load(dtmName);
else
    DTM = dtmName; %if matrix is given as input
end

%Pre-allocate
Data = -9999*ones(length(x)*length(y),3);

counter = 0;

for j = 1:length(y)
    for i = 1:length(x)
        counter = counter +1;
        x_temp = x(i);
        y_temp = y(j);
        
        Data(counter,1) = x_temp;
        Data(counter,2) = y_temp;
        
        d = find(DTM(:,1) == x_temp & DTM(:,2) == y_temp);
        
        % If  a point in the dtm file use the elevation else leave as -9999
        if ~isempty(d)
            Data(counter,3) = DTM(d,3);
        end
    end
end

% Write to a file
dlmwrite(outputName,Data,'delimiter',' ','precision',10);

%--------------------------------------------------------------------------