function [Gloms] = extract_all_gloms_SINGLE(xml_struct, img_ptr, verbose)
%extract_all_gloms_SINGLE extracts all gloms (in bounding boxes) within a
%single WSI, without image rescaling or downsampling
%   Inputs:
%       xml_struct:         Struct representing XML file
%       img_ptr:            Image pointer object from WSI library
%       verbose:            Logical, terminal output or not

% Pre-allocate Gloms struct
Gloms = cell(0);

% Traverse annotations
numAnnotations = length(xml_struct.Annotations.Annotation);

for Annotation_ticker = 1:numAnnotations
    Annotation = xml_struct.Annotations.Annotation(Annotation_ticker);
    annotationID = Annotation.Attributes.Id;
    
    % Traverse regions
    numRegions = length(Annotation.Regions.Region);
    if(verbose)
        fprintf('Reading Annotation %i/%i\n', Annotation_ticker, numAnnotations);
        fprintf('\tRegions found: %i\n', numRegions);
    end

    for Region_ticker = 1:numRegions
        Region = Annotation.Regions.Region{Region_ticker};
        regionID = Region.Attributes.Id;
        
        % Get vertex points, store in Nx2 matrix as [x,y] rows
        numVertices = length(Region.Vertices.Vertex);
        
        if numVertices < 3  % BAD RECORDING - SKIP THIS ONE
            warning('A%i.R%i SKIPPED - Less than three vertices.', Annotation_ticker, Region_ticker);
            continue
        end
        Vertices = nan(numVertices, 2);
        for Vertex_ticker = 1:numVertices
            Vertex = Region.Vertices.Vertex{Vertex_ticker};
            
            x_point = str2double(Vertex.Attributes.X);
            y_point = str2double(Vertex.Attributes.Y);
            Vertices(Vertex_ticker, 1:2) = [x_point, y_point];
        end
        
        % Convert vertices to polyshape, obtain bounding box and mask
        Shape = polyshape(Vertices);
        [BoundingBox_x, BoundingBox_y] = boundingbox(Shape);
        if isempty(BoundingBox_x)     % BAD RECORDING - SKIP THIS ONE
            warning('A%i.R%i SKIPPED - Empty bounding box.', Annotation_ticker, Region_ticker);
            continue
        end
        Vertices_norm = Vertices - [BoundingBox_x(1), BoundingBox_y(1)];
        Mask = poly2mask(Vertices_norm(:,1)', Vertices_norm(:,2)', ...
                              diff(BoundingBox_y), diff(BoundingBox_x));
        
        % Extract image of glom, and apply mask
        Image_Full = openslide_read_region(img_ptr, ...
            BoundingBox_x(1),    BoundingBox_y(1), ...
            diff(BoundingBox_x), diff(BoundingBox_y), ...
            0);
        Image_Full = Image_Full(:,:,2:4);
        
        % Save all data
        Gloms{length(Gloms)+1} = ...
                  struct('annotationID',  annotationID, ...
                         'regionID',      regionID, ...
                         'Vertices',      Vertices, ...
                         'BoundingBox',   [BoundingBox_x; BoundingBox_y], ...
                         'Vertices_norm', Vertices_norm, ...
                         'Mask',          Mask, ...
                         'Image_Full',    Image_Full);
    end     
    if(verbose)
        fprintf('\tRecording complete.\n');
    end
end

end

