% takes unthresholded corr matrix and returns thresholded binarized matrix,
% its degree distribution, and its entropy
function [w_bin_abs, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_abs)
   abs_thresh_matrix = threshold_absolute(unthresh_matrix, thresh_abs);
   w_bin_abs = weight_conversion(abs_thresh_matrix, 'binarize');
   ddist = degrees_und(w_bin_abs);
   
   entropy_bins = 120;
   f = figure('visible','off'); %initialize fig invisible
   hist1 = histogram(ddist,entropy_bins);
   hist1.Normalization = 'probability';
   binvalues = hist1.Values;
   nS = sum(-(binvalues(binvalues>0).*(log(binvalues(binvalues>0)))));
   clear f; clear hist1;
end