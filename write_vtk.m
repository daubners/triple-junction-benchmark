%% Writes matlab matrix as vtk output
%  Is limited to 2D, scalar data on a regular grid.
% Usage: write_vtk('example.vtk',field,0.1,0.1,'binary');

function write_vtk(filename,fieldname,field,dx,dy,format)
% Get the matrix dimensions.
[Nx,Ny,Nz] = size(field);

% Open the file.
out = fopen(filename, 'w');
if out == -1
    error('Cannot open file for writing.');
end

switch format
    case 'ascii'
        fprintf(out,'# vtk DataFile Version 3.0\n');
        fprintf(out,['Field variable ',fieldname,'\n']);
        fprintf(out,'ASCII\n');
        fprintf(out,'DATASET STRUCTURED_POINTS\n');
        fprintf(out,'DIMENSIONS %d %d %d\n',Nx,Ny,Nz);
        fprintf(out,'ORIGIN %d %d %d\n',-dx/2,-dy/2,0);
        fprintf(out,'SPACING %d %d %d\n',dx,dy,1);
        fprintf(out,'POINT_DATA %d\n',Nx*Ny*Nz);
        fprintf(out,['SCALARS ',fieldname,' float 1\n']);
        fprintf(out,'LOOKUP_TABLE default\n');
        fwrite(out, num2str(field(:)'));
    case 'binary'
        fprintf(out,'# vtk DataFile Version 3.0\n');
        fprintf(out,['Field variable ',fieldname,'\n']);
        fprintf(out,'BINARY\n');
        fprintf(out,'DATASET STRUCTURED_POINTS\n');
        fprintf(out,'DIMENSIONS %d %d %d\n',Nx,Ny,Nz);
        fprintf(out,'ORIGIN %d %d %d\n',-dx/2,-dy/2,0);
        fprintf(out,'SPACING %d %d %d\n',dx,dy,1);
        fprintf(out,'POINT_DATA %d\n',Nx*Ny*Nz);
        fprintf(out,['SCALARS ',fieldname,' float 1\n']);
        fprintf(out,'LOOKUP_TABLE default\n');
        fwrite(out, field(:),'float','ieee-be');
    otherwise
        error('Dataformat must be ascii or binary!');
end

% Close the file.
fclose(out);
