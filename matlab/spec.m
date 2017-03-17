function [P, nu, phases] = spec(tmseries, nu_min, nu_max, dnu, Pint)
% Calculate power density spectrum of timeseries.
%
% [P, nu, phases] = spec(tmseries)
% [P, nu, phases] = spec(tmseries, nu_min, nu_max, dnu, Pint)
%
% Written by
% Rasmus Handberg, rasmush@phys.au.dk
% Department of Physics and Astronomy, Aarhus University

	tmseries_isfile = ischar(tmseries);

	if tmseries_isfile
		if ~exist(tmseries, 'file')
			error('File not found');
		elseif nargin < 4
			data = importdata(tmseries);
			if isstruct(data)
				t = data.data(:,1);
			else
				t = data(:,1);
			end;
			clear data;
		end;
	elseif size(tmseries,2) ~= 2 && size(tmseries,2) ~= 3
		error('Invalid input');
	end;

	use_weights = ( size(tmseries,2)==3 && ~all(tmseries(:,3) == 1) );
	if nargin < 2
		nu_min = 0;
	end;
	if nargin < 3 || isempty(nu_max)
		nu_max = 0.5e6/median(diff(t));
	end;
	if nargin < 4 || isempty(dnu)
		dnu = 1e6/(max(t)-min(t));
	end;
	if nargin < 5 || isempty(Pint)
		Pint = dnu;
	end;

	% Save timeseries to file:
	if tmseries_isfile
		inputfile = tmseries;
	else
		A = tmseries.';
		inputfile = tempfile;
		save(inputfile, '-ascii', '-double', 'A');
		clear A;
	end;

	opts = '-quiet -tday -pd';
	
	% Save program input to file:
	cmd_inputfile = tempname;
	fid = fopen(cmd_inputfile, 'w');
	if use_weights
		fprintf(fid, '%d\n', use_weights);
	end;
	fprintf(fid, '%.16e\n', nu_min);
	fprintf(fid, '%.16e\n', nu_max);
	fprintf(fid, '%.16e\n', dnu);
	fprintf(fid, '%.16e\n', Pint);
	fclose(fid);

	% Run Fortran program:
	cmd = sprintf('!spec %s "%s" "%s" < "%s"', opts, inputfile, outputfile, cmd_inputfile);
	output = eval(cmd);

	% Delete the inputfile again:
	if ~tmseries_isfile, delete(inputfile); end;
	delete(cmd_inputfile);

	% Load the output:
	data = importdata(outputfile, ' ');
	if isstruct(data)
		nu = data.data(:,1);
		P = data.data(:,2);
		phases = data.data(:,3);
	else
		nu = data(:,1);
		P = data(:,2);
		phases = data(:,3);
	end;

	% Delete the outputfile again:
	delete(outputfile);
