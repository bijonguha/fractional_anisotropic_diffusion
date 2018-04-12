function w_thresh = wave_denoising(g, wlt, sorh,level)
x = uint8(g);
[thr,sorh,keepapp] = ddencmp('den','wv',x);
sorh = 's';
w_thresh = wdencmp('gbl',x,wlt,level,thr,sorh,keepapp);
end