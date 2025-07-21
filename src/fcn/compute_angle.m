function phi = compute_angle(v1, v2)

cos_phi = dot(v1, v2)/(norm(v1)*norm(v2));
% Sometimes truncation errors can make the cos value land slighly outside
% the domain causing an non real angle values, which in unacceptable.
if cos_phi < -1
    cos_phi = -1;
elseif cos_phi > 1
    cos_phi = 1;
end

phi = acos(cos_phi);


end