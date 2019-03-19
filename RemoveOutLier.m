function [Mat,Outlier] = RemoveOutLier(Mat)
out_lier_index = [];
Outlier =[];
for i = 1:length(Mat)
    Mat_t = Mat;
    Mat_t(i) = [];
    if (mean(Mat_t)/mean(Mat))<0.95
        out_lier_index = [out_lier_index , i];
    end
end

if ~isempty(out_lier_index)
    Outlier = Mat(out_lier_index);
    Mat(out_lier_index) = [];
end