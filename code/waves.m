function [zn] = waves(pr1, pr2, h, n)
    random = rand(1, n);
    i = 0;
    phi = random.*(2*pi);
    wstart = 0.3;
    wend = 1.4;
    dw = (wend - wstart)/n;
    for j = 1:n
        w(j) = wstart + (j - 1)*dw + dw/2;
        c(j) = quadl('spectrplotn', (wstart + (j - 1)*dw), (wstart + j*dw), 10e-5);
        c(j) = sqrt(2*c(j));
    end
    for t = pr1:h:pr2
        i = i + 1;
        ti(i) = t;
        zn(i) = 0;
        for j = 1:n
            zn(i) = zn(i) + (0.33)*c(j)*cos(w(j)*t + phi(1, j));
        end
    end
    %hold on
    %plot(ti, zn)
end
