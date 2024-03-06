function bicycle_tf= get_bicycle_no_fb_tf(h, V, a, b, c, lambda)
    g = 9.81; % graviational accelration
    bicycle_tf = tf([a*V*sin(lambda)/(b*h), (V^2*h-a*c*g)*sin(lambda)/(b*h^2)], [1, 0, -g/h]);
end