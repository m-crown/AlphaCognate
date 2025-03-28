import { useEffect, useRef } from "react";
import pako from "pako";

interface MolstarViewerProps {
  url: string;
}

const MolstarViewer: React.FC<MolstarViewerProps> = ({ url }) => {
  const viewerRef = useRef<HTMLDivElement | null>(null);
  console.log(url);
  useEffect(() => {
    if (!window.PDBeMolstarPlugin || !viewerRef.current || !url) return;

    let objectUrl: string | null = null; // Store generated object URL

    async function init() {
      try {
        // Fetch the gzipped data
        const response = await fetch(url);
        const compressedData = await response.arrayBuffer();

        // Decompress using pako
        const decompressedData = pako.ungzip(new Uint8Array(compressedData), { to: "string" });

        // Create a blob for Molstar
        const blob = new Blob([decompressedData], { type: "text/plain" });
        objectUrl = URL.createObjectURL(blob);

        // Initialize the Molstar plugin
        const viewerInstance = new PDBeMolstarPlugin();
        const options = {
          customData: { url: objectUrl, format: "cif" },
          bgColor: "black",
          hideControls: true,
          expanded: false,
          encoding: "cif",
        };

        viewerInstance.render(viewerRef.current, options);
      } catch (error) {
        console.error("Error loading gzip-compressed structure:", error);
      }
    }

    init();

    return () => {
      if (objectUrl) URL.revokeObjectURL(objectUrl); // Clean up object URL
      viewerRef.current?.replaceChildren(); // Cleanup viewer
    };
  }, [url]);

  return (
    <div
      ref={viewerRef}
      id="myViewer"
      style={{
        width: "100%",
        height: "100%",
        minWidth: "40vh",
        minHeight: "80vh",
        position: "relative",
      }}
    />
  );
};

export default MolstarViewer;
