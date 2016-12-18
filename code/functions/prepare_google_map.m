function IMG = prepare_google_map(X, Y, lonV, latV, img)

IMG        = zeros(size(X,1), size(X,2), 3);
IMG(:,:,1) = interp2(lonV, latV, img(:,:,1), X,Y);
IMG(:,:,2) = interp2(lonV, latV, img(:,:,2), X,Y);
IMG(:,:,3) = interp2(lonV, latV, img(:,:,3), X,Y);