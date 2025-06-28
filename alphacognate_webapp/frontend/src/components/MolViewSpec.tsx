import React, { useEffect, useRef } from "react";


interface MolViewSpecViewerProps {
    fileData: string | null;
  }

const MolViewSpecViewer: React.FC<MolViewSpecViewerProps> = ({ fileData }) => {
    if (!window.molstar || !fileData) return;
    const viewerRef = useRef<HTMLDivElement>(null);

    useEffect(() => {
        // Check that the Mol* global is loaded
        if (!window.molstar) {
        console.error("Mol* viewer not found");
        return;
        }

        const builder = window.molstar.PluginExtensions.mvs.MVSData.createBuilder();

        const structure = builder
        .download({ url: fileData })
        .parse({ format: "mmcif" })
        .modelStructure({})
        
        structure
            .component({ selector: "polymer" })
            .representation({ type: "cartoon" })
            .color({ color: "green" });

        structure
            .component({ selector: "ion" })
            .representation({ type: "ball_and_stick" })
            .color({custom: {molstar_color_theme_name: "element-index"}});

        const mvsData = builder.getState();
        
        const mvsj = window.molstar.PluginExtensions.mvs.MVSData.toMVSJ(mvsData);
        viewerRef.loadMvsData(mvsj, "mvsj", { replaceExisting: true });

        // Create the viewer
        window.molstar.Viewer.create(viewerRef.current!, {
        layoutIsExpanded: false,
        layoutShowControls: false,
        disabledExtensions: ['volseg'],
        }).then((viewer: any) => {
        viewer.loadMvsData(mvsj, "mvsj", { replaceExisting: true });
        });
    }, []);

    return (
        <div ref={viewerRef} style={{ 
            width: "100%", 
            height: "100%", 
            minWidth: "40vh", 
            minHeight: "80vh", 
            position: "relative" }} />
    );
    };

    export default MolViewSpecViewer;