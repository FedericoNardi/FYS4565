% SET PARTICLE NUMBER PARAMETER
function []= setParticleNumber( particleNumber )
    
    % set only particle number param
    params = getParams();
    params(1) = particleNumber;
    setParams(params);
    
end
