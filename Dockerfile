FROM mathworks/matlab:r2025b

# Switch to root to allow installing software into the system folders
USER root

# Download the MATLAB Package Manager (mpm) and install the toolbox
RUN wget -q https://www.mathworks.com/mpm/glnxa64/mpm && \
    chmod +x mpm && \
    ./mpm install \
        --release=R2025b \
        --destination=/opt/matlab/R2025b \
        --products Signal_Processing_Toolbox && \
    rm -f mpm /tmp/mathworks_root.log

# Switch back to the safe, default MATLAB user
USER matlab