function IMG = prepare_google_map(dem, lonV, latV, img)

IMG        = zeros(size(dem.X,1), size(dem.X,2), 3);
IMG(:,:,1) = interp2(lonV, latV, img(:,:,1), dem.X,dem.Y);
IMG(:,:,2) = interp2(lonV, latV, img(:,:,2), dem.X,dem.Y);
IMG(:,:,3) = interp2(lonV, latV, img(:,:,3), dem.X,dem.Y);