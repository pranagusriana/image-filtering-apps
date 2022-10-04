% Filter untuk notch filter pada citra uji A
filtersA = [
50 173 58 186;
65 198 78 206;
305 178 318 186;
325 198 333 211;
177 177 187 187;
197 197 207 207;
178 50 186 58;
198 70 206 78;
178 305 186 313;
198 325 206 333;
];

% Filter untuk notch filter pada citra uji B
filtersB = [
155 130 165 140;
090 130 100 140;
025 130 035 140;
355 130 365 140;
415 130 425 140;
485 130 495 140;
155 385 165 395;
090 385 100 395;
025 385 035 395;
355 385 365 395;
415 385 425 395;
485 385 495 395;
155 065 165 075;
090 065 100 075;
025 065 035 075;
355 065 365 075;
415 065 425 075;
485 065 495 075;
155 455 165 465;
090 455 100 465;
025 455 035 465;
355 455 365 465;
415 455 425 465;
485 455 495 465;
315 255 325 265;
370 255 380 265;
140 255 150 265;
200 255 210 265;
];

% Filter untuk notch filter pada citra uji C
filtersC = [
254 24 257 28;
236 32 238 52;
220 49 222 69;
204 66 206 86;
188 83 190 103;
172 100 174 120;
156 122 158 126;
124 154 126 158;
108 160 110 180;
92 177 94 197;
76 194 78 214;
60 211 62 231;
44 228 46 248;
25 252 7 255;
];

% Filter untuk notch filter pada citra uji D
filtersD = [
22 231 27 233;
22 331 27 333;
22 431 27 433;
35 231 43 233;
35 331 43 333;
35 431 43 433;
45 231 53 233;
45 331 53 333;
45 431 53 433;
92 231 102 233;
92 331 102 333;
92 431 102 433;
139 231 149 233;
139 331 149 333;
139 431 149 433;
169 231 179 233;
169 331 179 333;
169 431 179 433;
190 231 207 233;
190 331 207 333;
190 431 207 433;
218 231 229 233;
218 331 229 333;
218 431 229 433;
243 231 258 233;
243 431 258 433;
267 231 278 233;
267 331 278 333;
267 431 278 433;
289 231 305 233;
289 331 305 333;
289 431 305 433;
317 231 327 233;
317 331 327 333;
317 431 327 433;
347 231 357 233;
347 331 357 333;
347 431 357 433;
392 231 402 233;
392 331 402 333;
392 431 402 433;
442 231 453 233;
442 331 453 333;
442 431 453 433;
455 231 460 233;
455 331 460 333;
455 431 460 433;
467 231 472 233;
467 331 472 333;
467 431 472 433;
];


% FUngsi untuk melakukan notch filter
function out_img = notch_filter(Img, filters)
    [F, ~] = spatial2frequency(Img);
    F = fftshift(F);

    n = size(filters, 1);
    for i = 1:n
        [x1, y1, x2, y2] = deal(filters(i, 1), filters(i, 2), filters(i, 3), filters(i, 4));

        F(x1:x2, y1:y2) = 0;
    end

    out_img = frequency2spatial(ifftshift(F));
end

% Fungsi untuk melakukan highboost filter
function out_img = pass_highboost_filter(Img, filter_type, D0, n, A)
    [M, N, channels] = size(Img);
    %Step 1: Tentukan parameter padding, biasanya untuk citra f(x,y)
    % berukuran M x N, parameter padding P dan Q adalah P = 2M and Q = 2N.
    P = 2*M;
    Q = 2*N;
    %Step 2: Bentuklah citra padding fp(x,y) berukuran P X Q dengan
    % menambahkan pixel-pixel bernilai nol pada f(x, y).
    % Proses padding dilakukan dengan menambahkan sebanyak M elemen 0 pada rows dan N
    % elemen 0 cols
    fp = padarray(Img, [M, N], 'post');
    %Step 3: Bangkitkan fungsi penapis H berukuran P x Q
    if strcmp(filter_type, 'hbilpf')
        H = A - ILPF(D0, P, Q);
    elseif strcmp(filter_type, 'hbglpf')
        H = A - GLPF(D0, P, Q);
    elseif strcmp(filter_type, 'hbblpf')
        H = A - BLPF(D0, n, P, Q);
    end
    % Untuk setiap channels warna
    out_img = zeros(M, N, channels);
    for c = 1:channels
        %Step 4: Lakukan transformasi Fourier pada fpad(x, y)
        [F, ~] = spatial2frequency(fp(:, :, c));
        %Step 5: Kalikan F dengan H
        G = H.*F;
        %Step 6: Ambil bagian real dari inverse FFT of G:
        g = frequency2spatial(G); % apply the inverse, discrete Fourier transform
        %Step 7: Potong bagian kiri atas sehingga menjadi berukuran citra semula
        out_img(:, :, c) = g(1:M, 1:N); % Resize the image to undo padding
    end
    out_img = uint8(out_img);
end

% Create mean filter for spatial domain with nxn size
function filter = create_mean_filter(n)
    filter = ones(n);
    filter = filter/sum(filter, "all");
end

% create gaussian filter for spatial domain with nxn size
function filter = create_gaussian_filter(n, sigma)
    %// Generate horizontal and vertical co-ordinates, where
    %// the origin is in the middle
    ind = -floor(n/2) : floor(n/2);
    [X, Y] = meshgrid(ind, ind);
    
    %// Create Gaussian Mask
    h = (1 / (2 * pi * sigma ^ 2)) * exp(-(X.^2 + Y.^2) / (2*sigma*sigma));
    
    %// Normalize so that total area (sum of all weights) is 1
    filter = h / sum(h(:));
end

% Fungsi untuk melakukan passing low dan high pass filter
function out_img = pass_filter(Img, filter_type, D0, n)
    [M, N, channels] = size(Img);
    %Step 1: Tentukan parameter padding, biasanya untuk citra f(x,y)
    % berukuran M x N, parameter padding P dan Q adalah P = 2M and Q = 2N.
    P = 2*M;
    Q = 2*N;
    %Step 2: Bentuklah citra padding fp(x,y) berukuran P X Q dengan
    % menambahkan pixel-pixel bernilai nol pada f(x, y).
    % Proses padding dilakukan dengan menambahkan sebanyak M elemen 0 pada rows dan N
    % elemen 0 cols
    fp = padarray(Img, [M, N], 'post');
    %Step 3: Bangkitkan fungsi penapis H berukuran P x Q
    if strcmp(filter_type, 'ilpf')
        H = ILPF(D0, P, Q);
    elseif strcmp(filter_type, 'glpf')
        H = GLPF(D0, P, Q);
    elseif strcmp(filter_type, 'blpf')
        H = BLPF(D0, n, P, Q);
    elseif strcmp(filter_type, 'ihpf')
        H = IHPF(D0, P, Q);
    elseif strcmp(filter_type, 'ghpf')
        H = GHPF(D0, P, Q);
    elseif strcmp(filter_type, 'bhpf')
        H = BHPF(D0, n, P, Q);
    end
    % Untuk setiap channels warna
    out_img = zeros(M, N, channels);
    for c = 1:channels
        %Step 4: Lakukan transformasi Fourier pada fpad(x, y)
        [F, ~] = spatial2frequency(fp(:, :, c));
        %Step 5: Kalikan F dengan H
        G = H.*F;
        %Step 6: Ambil bagian real dari inverse FFT of G:
        g = frequency2spatial(G); % apply the inverse, discrete Fourier transform
        %Step 7: Potong bagian kiri atas sehingga menjadi berukuran citra semula
        out_img(:, :, c) = g(1:M, 1:N); % Resize the image to undo padding
    end
    out_img = uint8(out_img);
end

% Fungsi untuk mendapatkan kernel BHPF (Butterworth High Pass Filter
function H = BHPF(D0, n, P, Q)
    H = 1 - BLPF(D0, n, P, Q);
end

% Fungsi untuk mendapatkan kernel GHPF (Gaussian High Pass Filter)
function H = GHPF(D0, P, Q)
    H = 1 - GLPF(D0, P, Q);
end

% Fungsi untuk mendapatkan kernel IHPF (Ideal High Pass Filter)
function H = IHPF(D0, P, Q)
    H = 1 - ILPF(D0, P, Q);
end

% Fungsi untuk mendapatkan kernel BLPF (Butterworth Low Pass Filter)
function H = BLPF(D0, n, P, Q)
    D = designing_filter(P, Q);
    H = 1./(1 + (D./D0).^(2*n));
end

% Fungsi untuk mendapatkan kernel GLPF (Gaussian Low Pass Filter)
function H = GLPF(D0, P, Q)
    D = designing_filter(P, Q);
    H = exp(-(D.^2)./(2*(D0^2)));
end

% Fungsi untuk mendapatkan kernel ILPF (Ideal Low Pass Filter)
function H = ILPF(D0, P, Q)
    D = designing_filter(P, Q);
    H = double(D <= D0);
end

% Fungsi template untuk membuat filter pada ranah frequency
function D = designing_filter(P, Q)
    % Set up range of variables.
    u = 0:(P-1);
    v = 0:(Q-1);
    % Compute the indices for use in meshgrid
    idx = find(u > P/2);
    u(idx) = u(idx) - P;
    idy = find(v > Q/2);
    v(idy) = v(idy) - Q;
    % Compute the meshgrid arrays
    % MATLAB library function meshgrid(v, u) returns
    % 2D grid which contains the coordinates of vectors
    % v and u. Matrix V with each row is a copy 
    % of v, and matrix U with each column is a copy of u
    [V, U] = meshgrid(v, u);
    D = sqrt(U.^2 + V.^2);
end

% Fungsi untuk mengubah spatial domain ke frequency domain
function [ft_img, fs_img] = spatial2frequency(Img)
    ft_img = fft2(double(Img));
    fs_img = uint8(log(1+abs(fftshift(ft_img))));
end

% Fungsi untuk mengubah frequency domain ke spatial domain
function Img = frequency2spatial(ft_img)
    Img = uint8(real(ifft2(ft_img)));
end

% Fungsi untuk melakukan padding pada image agar output shape hasil
% konvolusi sama dengan input shape image
function Img = padsameimg(Img, kernel_height, kernel_width)
    two_padsize_h = kernel_height - 1;
    two_padsize_w = kernel_width - 1;
    if mod(two_padsize_h, 2) == 1 && mod(two_padsize_w, 2) == 1
        padsize_h_1 = floor(two_padsize_h/2);
        padsize_h_2 = floor(two_padsize_h/2) + 1;
        padsize_w_1 = floor(two_padsize_w/2);
        padsize_w_2 = floor(two_padsize_w/2) + 1;
        Img = double(padarray(Img, [padsize_h_1, padsize_w_1], 0, 'pre'));
        Img = double(padarray(Img, [padsize_h_2, padsize_w_2], 0, 'post'));
    elseif mod(two_padsize_h, 2) == 1 && mod(two_padsize_w, 2) == 0
        padsize_h_1 = floor(two_padsize_h/2);
        padsize_h_2 = floor(two_padsize_h/2) + 1;
        padsize_w_1 = two_padsize_w/2;
        padsize_w_2 = two_padsize_w/2;
        Img = double(padarray(Img, [padsize_h_1, padsize_w_1], 0, 'pre'));
        Img = double(padarray(Img, [padsize_h_2, padsize_w_2], 0, 'post'));
    elseif mod(two_padsize_h, 2) == 0 && mod(two_padsize_w, 2) == 1
        padsize_h_1 = two_padsize_h/2;
        padsize_h_2 = two_padsize_h/2;
        padsize_w_1 = floor(two_padsize_w/2);
        padsize_w_2 = floor(two_padsize_w/2) + 1;
        Img = double(padarray(Img, [padsize_h_1, padsize_w_1], 0, 'pre'));
        Img = double(padarray(Img, [padsize_h_2, padsize_w_2], 0, 'post'));
    elseif mod(two_padsize_h, 2) == 0 && mod(two_padsize_w, 2) == 0
        padsize_h_1 = two_padsize_h/2;
        padsize_h_2 = two_padsize_h/2;
        padsize_w_1 = two_padsize_w/2;
        padsize_w_2 = two_padsize_w/2;
        Img = double(padarray(Img, [padsize_h_1, padsize_w_1], 0, 'pre'));
        Img = double(padarray(Img, [padsize_h_2, padsize_w_2], 0, 'post'));
    end
end

% Same padding, output image shape = input image shape
function conv_img = conv2d(Img, kernel)
    Img = double(Img); % Samakan format image yang diterima
    [kernel_height, kernel_width] = size(kernel); % Dapatkan size dari kernel
    Img = padsameimg(Img, kernel_height, kernel_width); % Lakukan padding pada image agar output shape hasil konvolusi sama dengan input shapenya
    [height, width, channels] = size(Img); % Dapatkan size dari image
    [output_height, output_width] = deal(height - kernel_height + 1, width - kernel_width + 1); % Calculate output shape
    conv_img = zeros(output_height, output_width, channels); % Initiate image output hasil convolution
    % Untuk setiap channel warna gambar lakukan konvolui
    for c = 1:channels
        conv_img(:, :, c) = conv_operation(Img(:, :, c), kernel);
    end
    % Ubah ke bentuk uint8
    conv_img = uint8(conv_img);
end

function conv_img = conv_operation(mat, kernel)
    [height, width] = size(mat);
    [kernel_height, kernel_width] = size(kernel);
    [output_height, output_width] = deal(height - kernel_height + 1, width - kernel_width + 1);
    conv_img = zeros(output_height, output_width);
    % Operasi konvolusi
    for i = 1:output_height
        for j = 1:output_width
            curwindow = mat(i:i+kernel_height-1, j:j+kernel_width-1);
            conv_img(i, j) = sum(curwindow .* kernel, "all");
        end
    end
end

% Same padding, output image shape = input image shape
function medfil_img = medfil(Img, window_size)
    Img = double(Img); % Samakan format image yang diterima
    [window_height, window_width] = deal(window_size(1), window_size(2)); % Dapatkan window size
    Img = padsameimg(Img, window_height, window_width); % Lakukan padding pada image agar output shape hasil medfilter sama dengan input shapenya
    [height, width, channels] = size(Img); % Dapatkan size image
    [output_height, output_width] = deal(height - window_height + 1, width - window_width + 1); % Calculate output shape
    medfil_img = zeros(output_height, output_width, channels); % Initiate image output hasil medfilter
    % Untuk setiap channel warna gambar
    for c = 1:channels
        medfil_img(:, :, c) = medfil_operation(Img(:, :, c), window_size); % Lakukan operasi medfilter
    end
    % Ubah ke bentuk uint8
    medfil_img = uint8(medfil_img);
end

function medfil_img = medfil_operation(mat, window_size)
    [height, width] = size(mat);
    [window_height, window_width] = deal(window_size(1), window_size(2));
    [output_height, output_width] = deal(height - window_height + 1, width - window_width + 1);
    medfil_img = zeros(output_height, output_width);
    for i = 1:output_height
        for j = 1:output_width
            curwindow = mat(i:i+window_height-1, j:j+window_width-1);
            medfil_img(i, j) = median(curwindow, "all"); % Operasi Median Filter
        end
    end
end