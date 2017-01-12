
choiceGrid='v4';

%need to replace the missing v3 files...
%
%choose grid for test
% fprintf('\n Please set the demo grid : \n');
% fprintf('    ''llc'': lat-lon-cap grid (5 faces). \n');
% fprintf('    ''ll'' : simple lat-lon grid (1 face). \n');
% choiceGrid=input(' and type return. ''llc'' is the default.\n');
% if isempty(choiceGrid); choiceGrid='llc'; end;
% if strcmp(choiceGrid,'llc'); choiceGrid='v4';
% elseif strcmp(choiceGrid,'ll'); choiceGrid='v3';
% else; error('wrong grid choice');
% end;

%choose verbise level
fprintf('\n Please set the amount of explanatory text display :\n');
fprintf('    0: none.\n');
fprintf('    1: comments.\n');
fprintf('    2: comments preceeded with calling sequence.\n');
fprintf('    3: same as 2, but preceeded with pause.\n');
verbose=input(' and/or type return. 0 is the default. \n');
if isempty(verbose); verbose=0; end;

%initialize environment variables and mygrid
gcmfaces_global;
myenv.verbose=verbose;
if myenv.verbose>0;
    gcmfaces_msg('* set path and environment variables (myenv) by calling gcmfaces_global');
end;

fprintf('\n\n');
gcmfaces_msg('/////////////////////////////////////////');
gcmfaces_msg('demo of transport computations','// PART 1 :');
gcmfaces_msg('please hit return','// >> ');
gcmfaces_msg('/////////////////////////////////////////');
fprintf('\n\n');
pause;
diags=example_transports(choiceGrid);
if ~isempty(diags); 
    example_transports_disp(diags);
end;

fprintf('\n\n');
gcmfaces_msg('/////////////////////////////////////');
gcmfaces_msg('demo of plotting routines','// PART 2 :');
gcmfaces_msg('please hit return','// >> ');
gcmfaces_msg('/////////////////////////////////////');
fprintf('\n\n');
pause;
plot_one_field(choiceGrid,0);%%incl. call to m_map_gcmfaces;

fprintf('\n\n');
gcmfaces_msg('///////////////////////////////////////////');
gcmfaces_msg('demo of gridding v interpolating','// PART 3 :');
gcmfaces_msg('please hit return','// >> ');
gcmfaces_msg('///////////////////////////////////////////');
fprintf('\n\n');
pause;
example_bin_average(choiceGrid,1);%incl. call to example_smooth.m;
example_griddata(choiceGrid);
example_interp(choiceGrid);
example_faces2latlon2faces(choiceGrid);



