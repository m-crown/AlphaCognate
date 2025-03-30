import React from "react";

interface FocusButtonProps {
  viewerInstanceRef: React.MutableRefObject<any>;
}

const FocusButton: React.FC<FocusButtonProps> = ({ viewerInstanceRef }) => {
  const handleClick = () => {
    console.log(viewerInstanceRef)
    if (viewerInstanceRef.current) {
      viewerInstanceRef.current.visual.select({
        data: [{ struct_asym_id: "M", color: { r: 255, g: 255, b: 0 }, focus: true }],
      });
    }
  };
  const resetClick = () => {
    if (viewerInstanceRef.current) {
        viewerInstanceRef.current.visual.reset({highlightColor: true })
    }};

  return <button onClick={handleClick} onDoubleClick={resetClick}>Select & Focus on Chain A</button>;
};

export default FocusButton;