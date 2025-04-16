import { useEffect, useRef } from "react";
import pako from "pako";

interface MolstarViewerProps {
  url: string;
  viewerInstanceRef: React.RefObject<any>;
}

const MolstarViewer: React.FC<MolstarViewerProps> = ({ url, viewerInstanceRef }) => {
  const viewerRef = useRef<HTMLDivElement | null>(null);

  useEffect(() => {
    if (!window.PDBeMolstarPlugin || !viewerRef.current || !url) return;

    let objectUrl: string | null = null;

    async function init() {
      try {
        const response = await fetch(url);
        const compressedData = await response.arrayBuffer();
        const decompressedData = pako.ungzip(new Uint8Array(compressedData), { to: "string" });
        objectUrl = URL.createObjectURL(new Blob([decompressedData], { type: "text/plain" }));

        viewerInstanceRef.current = new window.PDBeMolstarPlugin();
        viewerInstanceRef.current.render(viewerRef.current, {
          customData: { url: objectUrl, format: "cif" },
          bgColor: "black",
          hideControls: true,
          expanded: false,
          encoding: "cif",
        });
      } catch (error) {
        console.error("Error loading structure:", error);
      }
    }

    init();

    return () => {
      if (objectUrl) URL.revokeObjectURL(objectUrl);
      viewerRef.current?.replaceChildren();
    };
  }, [url]);

  return <div ref={viewerRef} style={{ width: "100%", height: "100%", minWidth: "40vh", minHeight: "80vh", position: "relative" }} />;
};

export default MolstarViewer;
